"""
lumip_utils.py
================
Shared helper functions for the `code/revision` LUMIP afforestation/reforestation
analysis notebooks (ssp370 vs ssp126Lu, LUMIP deforestation/afforestation
experiment).

Every notebook in `code/revision` that needs one of these functions should
import it from here (`from lumip_utils import ...`) instead of carrying its own
copy-pasted definition. This is the single source of truth for two corrections
made during the August 2026 revision:

1. **Evaporative-regime classification (`regime_limit` / `classify_regime`)**
   was computed on the *raw* monthly hfls/pr/ts series. Because all three fields
   have a strong shared seasonal cycle, the water-vs-energy-limited correlation
   test was almost entirely driven by that shared seasonality rather than by the
   interannual land-atmosphere coupling signal it is meant to capture (the
   Seneviratne/Denissen definition of the evaporative regime is based on
   interannual anomalies). A deseasonalizing fix (removing the monthly
   climatology, `deseason`, from hfls/pr/ts before correlating) was added and is
   still available as `method="deseason"` below -- **but by the user's explicit
   request (Aug 2026) the default has been reverted to `method="raw"`**, i.e.
   the original, non-deseasonalized behaviour, everywhere in `code/revision`.
   `regime_limit`/`classify_regime`'s default argument controls this for every
   notebook that calls them without passing `method=`; the two notebooks with
   their own explicit `REGIME_METHOD` config cell (`water_avail_bootstrap.ipynb`,
   `water_avail_final.ipynb`) have that variable set to `"raw"` for the same
   reason. See `REVISION_NOTES.md` for the full history of this back-and-forth.

2. **4GST phenology classification (`compute_4gst`, replacing the old
   `classify_phenology` / `phen_4gst`)** had a bug in the evergreen test: it
   compared the seasonal amplitude to the annual *minimum* LAI
   (`amp < 0.25 * lai_min`) instead of the annual *mean* LAI
   (`amp / lai_mean < 0.25`), and had no guard against `lai_mean <= 0`. Since
   `lai_min` is frequently close to zero (any cell with a bare-soil month),
   `0.25 * lai_min` is itself close to zero and the evergreen branch almost
   never triggered, so a nontrivial fraction of weakly-seasonal/evergreen
   cells were being run through the onset/offset slope classifier instead
   and got a spurious onset/offset window instead of "all 12 months are
   growing season". `compute_4gst` below is the corrected version, supplied
   directly by the user (Peano et al. 2019 4GST method); it is used everywhere
   growing-season classification is needed and its outputs (onset, offset,
   gst_type) are cached to disk by `grow_seas.ipynb` for the other notebooks
   to load (see `growing_season_mask` below and `grow_seas.ipynb`).

Everything else in this module (`lag_linregress_3D`, `cell_weight`, `lon180`,
`xr_fix_time`, `xr_clean`, `detrend_dim`, `multi_apply_along_axis`,
`xr_multipletest`, `xr_multimodel_sign`) is copied verbatim from the original
notebooks (byte-for-byte identical logic) purely so every revised notebook can
`import` a single copy instead of six slightly-drifted ones.
"""
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import linregress


# ---------------------------------------------------------------------------
# Generic grid / regridding / multi-model helpers (unchanged from the
# original notebooks -- copied here only to remove duplication)
# ---------------------------------------------------------------------------

def cell_weight(ds):
    """Grid-cell area (m^2) on a regular lat-lon grid, up to a constant that
    cancels in any weighted mean/sum. Identical across all original notebooks."""
    R = 6.371e6
    dphi = np.deg2rad(ds.lat[1] - ds.lat[0])
    dlam = np.deg2rad(ds.lon[1] - ds.lon[0])
    dlat = R * dphi * xr.ones_like(ds.lon)
    dlon = R * dlam * np.cos(np.deg2rad(ds.lat))
    return dlon * dlat


def lon180(ds):
    """Convert a 0-360 longitude coordinate to -180-180 and sort by it."""
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds.sortby(ds.lon)


def xr_fix_time(xr_in, date_start, date_end):
    """Replace a (possibly non-standard-calendar) time coordinate with a
    regular monthly gregorian range, so different ESM calendars align."""
    if xr_in["time"].dt.calendar in ("noleap", "360_day"):
        xr_in = xr_in.convert_calendar(calendar="gregorian", align_on="date")
    time = pd.date_range(start=date_start, end=date_end, freq="ME").to_numpy(dtype="datetime64[ns]")
    return xr_in.assign_coords(time=time)


def xr_clean(xr_in, dims):
    """Drop bookkeeping coordinates/dims (bnds, member_id, ...) that differ
    across ESMs and are not needed for the analysis."""
    data = xr_in.copy()
    for d in dims:
        if d in data.coords and d in data.dims:
            data = data.drop_vars(d)
            data = data.isel({d: 0})
        if d in data.coords:
            data = data.drop_vars(d)
        if d in data.dims and d not in data.coords:
            data = data.isel({d: 0})
        if data.attrs != {} and d in data.data_vars:
            data = data.drop_vars(d)
    return data


def detrend_dim(da, dim, degree):
    """Remove a polynomial trend of the given degree along `dim`."""
    p = da.polyfit(dim=dim, deg=degree)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit


def lag_linregress_3D(x, y, lagx=0, lagy=0):
    """
    Grid-cell OLS of y on x along the (aligned) `time` dimension, with optional
    lags. x, y: xr.DataArrays of any shape with `time` as the first dimension.
    Returns a dict with "cor" (Pearson correlation), "coef" (regression slope)
    and "pval" (two-sided p-value of the slope/correlation).
    """
    x, y = xr.align(x, y)

    if lagx != 0:
        x = x.shift(time=-lagx).dropna(dim="time")
        x, y = xr.align(x, y)
    if lagy != 0:
        y = y.shift(time=-lagy).dropna(dim="time")
        x, y = xr.align(x, y)

    n = x.shape[0]
    xmean, ymean = x.mean(axis=0), y.mean(axis=0)
    xstd, ystd = x.std(axis=0), y.std(axis=0)

    cov = np.sum((x - xmean) * (y - ymean), axis=0) / n
    cor = cov / (xstd * ystd)
    slope = cov / (xstd ** 2)

    tstats = cor * np.sqrt(n - 2) / np.sqrt(1 - cor ** 2)
    from scipy.stats import t
    pval = t.sf(tstats, n - 2) * 2
    pval = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

    return {"cor": cor, "coef": slope, "pval": pval}


def multi_apply_along_axis(func1d, axis, arrs, *args, **kwargs):
    """Apply `func1d(A, B, C, ...)` along `axis` to N-D arrays in `arrs`."""
    carrs = np.concatenate(arrs, axis)
    offsets = []
    start = 0
    for i in range(len(arrs) - 1):
        start += arrs[i].shape[axis]
        offsets.append(start)

    def helperfunc(a, *args, **kwargs):
        split = np.split(a, offsets)
        return func1d(*[*split, *args], **kwargs)

    return np.apply_along_axis(helperfunc, axis, carrs, *args, **kwargs)


def xr_multipletest(p, alpha=0.05, method="fdr_bh", **multipletests_kwargs):
    """Apply statsmodels.stats.multitest.multipletests across an xr.object."""
    from statsmodels.stats.multitest import multipletests
    p_stacked = p.stack(s=p.dims)
    mask = np.isfinite(p_stacked)
    pvals_corrected = np.full(p_stacked.shape, np.nan)
    reject = np.full(p_stacked.shape, np.nan)
    reject[mask] = multipletests(p_stacked[mask], alpha=alpha, method=method, **multipletests_kwargs)[0]
    pvals_corrected[mask] = multipletests(p_stacked[mask], alpha=alpha, method=method, **multipletests_kwargs)[1]

    def unstack(vals, template):
        out = template.copy()
        out.values = vals
        return out.unstack()

    return unstack(reject, p_stacked), unstack(pvals_corrected, p_stacked)


def xr_multimodel_sign(xr_in, models):
    """Multi-model mean, std, and sign-agreement fraction across the 'esm' dim."""
    mm = xr.concat(xr_in, dim="esm", coords="minimal", compat="override")
    std = mm.std(dim="esm")
    sign = xr.where(mm > 0, 1, -1)
    agreement = np.abs(sign.sum(dim="esm") / len(models))
    mean = mm.mean(dim="esm")
    return mean, std, agreement


# ---------------------------------------------------------------------------
# FIX 1 -- evaporative-regime classification (deseasonalized)
# ---------------------------------------------------------------------------

def deseason(da):
    """Remove the monthly climatology, leaving interannual anomalies.
    `da.groupby("time.month") - da.groupby("time.month").mean("time")`, tidied
    of the leftover `month` coordinate and returned with a canonical dim order.
    """
    out = da.groupby("time.month") - da.groupby("time.month").mean("time")
    return out.drop_vars("month", errors="ignore").transpose("time", ...)


def classify_regime(ds, method="raw"):
    """1 = water-limited (LE covaries more with Pr), 0 = energy-limited (more
    with Ts), NaN = undefined/ocean. 
    `method="raw"` default correlates the raw monthly hfls/pr/ts series; 
    `method="deseason"` removes the shared seasonal cycle first;

    NaN is preserved explicitly (unlike the original `xr.where(cor_pr > cor_ts,
    1, 0)`, which silently mapped every undefined/ocean cell to energy-limited).
    """
    le, pr, ts = ds["hfls"], ds["pr"], ds["ts"]
    if method == "deseason":
        le, pr, ts = deseason(le), deseason(pr), deseason(ts)
    cor_pr = lag_linregress_3D(le, pr)["cor"].compute()
    cor_ts = lag_linregress_3D(le, ts)["cor"].compute()
    ok = cor_pr.notnull() & cor_ts.notnull()
    return xr.where(ok, xr.where(cor_pr > cor_ts, 1.0, 0.0), np.nan)


def regime_limit(dataset, method="raw", drop=False):
    """Same call signature and return contract as the original
    `regime_limit(dataset)` used throughout the pre-revision notebooks (a pair
    of 2-D `lai` DataArrays, valid only where water-/energy-limited
    respectively), so downstream code (`mask_aff_water[m].notnull()`, etc.) is
    unaffected regardless of `method`.

    `method="raw"` (default) correlates the raw monthly hfls/pr/ts series --
    this reproduces the original, pre-revision behaviour exactly, restored as
    the default by explicit user request (Aug 2026) after briefly defaulting to
    `method="deseason"` (still available, and statistically the more defensible
    choice -- it removes the shared seasonal cycle from hfls/pr/ts before
    correlating, see module docstring -- but not the current default here).

    `drop` matches the original per-notebook choice: `grow_seas`/`water_avail*`/
    `seb_grow_seas`/`roughness_scale*` used `drop=False` (keep full grid, NaN
    outside the regime); `mlr_comp`/`mlr_proc` used `drop=True` (compact away
    all-NaN slices). Both are label-aligned by xarray downstream regardless of
    `drop`, but each notebook here keeps its original setting.
    """
    le, pr, ts = dataset["hfls"], dataset["pr"], dataset["ts"]
    if method == "deseason":
        le, pr, ts = deseason(le), deseason(pr), deseason(ts)

    cor_pr = lag_linregress_3D(le, pr)["cor"].compute()
    cor_ts = lag_linregress_3D(le, ts)["cor"].compute()

    water_limited = dataset.isel(time=-1).where(cor_pr > cor_ts, drop=drop)["lai"]
    energy_limited = dataset.isel(time=-1).where(cor_pr < cor_ts, drop=drop)["lai"]
    return water_limited, energy_limited


# ---------------------------------------------------------------------------
# FIX 2 -- 4GST growing-season-type classification (evergreen bug fixed)
# ---------------------------------------------------------------------------

def compute_4gst(clim_da, month_dim="month"):
    """
    Growing-season-type classification (4GST method, Peano et al. 2019),
    applied per-pixel to the climatological (12-month) mean LAI cycle.
    Reference implementation as supplied directly by the user; GST values:
      1, 2 -- single growing season (unimodal), mirror-image slope patterns
      3    -- two growing seasons per year (double-cropping / bimodal)
      4    -- evergreen / weakly seasonal (amplitude below 25% of the annual
              mean LAI); all 12 months are retained as growing season
      0    -- undefined / fallback

    Corrects the evergreen test relative to the original `classify_phenology`
    / `phen_4gst` found across the pre-revision notebooks, which tested
    `amp < 0.25 * lai_min` (amplitude vs. the annual *minimum*, with no floor
    on `lai_min`) instead of `amp / lai_mean < 0.25` (amplitude vs. the annual
    *mean*, guarded against `lai_mean <= 0`). Because `lai_min` is close to
    zero for most vegetated cells, the old test's threshold collapsed toward
    zero and rarely fired, so evergreen/weakly-seasonal cells were frequently
    misrouted into the onset/offset slope classifier instead of being flagged
    as evergreen (all-year growing season).
    """
    def _pixel(series):
        series = np.asarray(series)
        if np.all(np.isnan(series)):
            return (np.nan, np.nan, np.nan)
        midpos = np.nanargmax(series)
        rotated = np.roll(series, 6 - midpos)
        lai_min, lai_max = np.nanmin(series), np.nanmax(series)
        lai_mean = np.nanmean(series)
        amp = lai_max - lai_min
        above = np.where(series > lai_min + 0.2 * amp)[0]
        # evergreen: seasonal amplitude small relative to the mean annual LAI
        if lai_mean > 0 and amp / lai_mean < 0.25:
            return (1, 12, 4)
        if above.size == 0:
            return (1, 12, 0)
        slopes = []
        for off in (0, 3, 6, 9):
            seg = np.nan_to_num(rotated[off:off + 3], nan=0.0)
            slopes.append(linregress(np.arange(3), seg).slope)
        mapping = {(1, 1, -1, -1): 1, (-1, -1, 1, 1): 2, (1, -1, 1, -1): 3, (-1, 1, -1, 1): 3}
        return (above[0] + 1, above[-1] + 1,
                mapping.get(tuple(np.sign(slopes).astype(int)), 0))

    def _chunk(arr):
        # arr: (n_lat, n_lon, 12)  [core dim last via apply_ufunc]
        n_lat, n_lon, _ = arr.shape
        start = np.full((n_lat, n_lon), np.nan)
        end = np.full((n_lat, n_lon), np.nan)
        gst = np.full((n_lat, n_lon), np.nan)
        for i in range(n_lat):
            for j in range(n_lon):
                start[i, j], end[i, j], gst[i, j] = _pixel(arr[i, j, :])
        return start, end, gst

    start_da, end_da, gst_da = xr.apply_ufunc(
        _chunk, clim_da,
        input_core_dims=[[month_dim]],
        output_core_dims=[[], [], []],
        vectorize=False,
        dask="parallelized",
        output_dtypes=[np.float64, np.float64, np.float64],
    )
    return xr.Dataset({"gst_start_month": start_da, "gst_end_month": end_da, "gst_type": gst_da})


def growing_season_reduce(series, months, onset, offset):
    """
    Corrected drop-in replacement for the original `growing_season()` /
    `phen_gs()` scalar reducers found in `grow_seas.ipynb`: the `nanmean` of
    `series` restricted to the months between `onset` and `offset`
    (inclusive), handling the Southern-Hemisphere wrap-around case
    (onset > offset) explicitly.

    The original mask, `(months >= onset) & (months <= offset)`, silently
    evaluates to all-False for *every* month whenever onset > offset (a
    growing season that crosses the calendar-year boundary, e.g.
    onset=October, offset=March) -- i.e. it masked out the entire year,
    including the actual growing season, at every such pixel. This is the
    "must mask only the portions of the year where vegetation is active" bug:
    for wrap-around pixels the old code masked none of the active months (and
    all of the inactive ones), so `np.nanmean(series[mask])` silently reduced
    to a mean of an all-empty selection.

    Meant for use with `xr.apply_ufunc(growing_season_reduce, series, months,
    onset, offset, input_core_dims=[["time"], ["time"], [], []],
    output_core_dims=[[]], vectorize=True, ...)`, exactly as the original
    `growing_season()` was used -- i.e. this collapses the selected period to
    a single climatological growing-season-mean value per pixel. For an
    analysis that instead needs the growing season masked *without* collapsing
    the time dimension (e.g. an interannual regression), use
    `growing_season_mask` below with `.where(...)` instead.
    """
    if np.isnan(onset) or np.isnan(offset):
        return np.nanmean(series)
    if onset <= offset:
        mask = (months >= onset) & (months <= offset)
    else:
        mask = (months >= onset) | (months <= offset)
    if not np.any(mask):
        return np.nan
    return np.nanmean(series[mask])


def growing_season_mask(time_da, onset, offset):
    """Boolean (time, lat, lon) mask, True inside the growing season defined by
    per-pixel `onset`/`offset` month numbers (1-12), False outside it and at
    pixels with undefined onset/offset. Handles the southern-hemisphere
    wrap-around case (onset > offset) explicitly.

    `time_da` is any DataArray/coordinate exposing `.dt.month` along `time`
    (e.g. `ds.time`). Unlike the pre-revision `growing_season()` /
    `phen_gs()` helpers -- which collapsed the whole selected period to a
    single `nanmean` scalar per pixel, discarding the time dimension -- this
    returns a mask the same shape as the data, so that
    `ds.where(growing_season_mask(ds.time, onset, offset))` keeps the full
    (annual/monthly) time resolution and simply sets non-growing-season
    time steps to NaN. This is what "masking only the portions of the year
    where vegetation is active" means: downstream notebooks choose their own
    reduction (annual mean/sum, climatology, regression, ...) on the masked,
    still time-resolved data, rather than getting a pre-collapsed scalar.
    """
    months = time_da.dt.month
    cond_normal = (months >= onset) & (months <= offset)
    cond_crossed = (months >= onset) | (months <= offset)
    mask = xr.where(onset <= offset, cond_normal, cond_crossed)
    return mask.where(np.isfinite(onset) & np.isfinite(offset), False)


# ---------------------------------------------------------------------------
# Growing-season cache I/O contract shared by grow_seas.ipynb (writer) and
# every other revised notebook (readers) -- see grow_seas.ipynb for how these
# are produced.
# ---------------------------------------------------------------------------

GS_CACHE_DIR = "./cache/growing_season"


def gs_cache_path(scenario, model, cache_dir=GS_CACHE_DIR):
    """scenario: 'ctl' or 'aff'. Path to the cached phenology (onset, offset,
    gst_type) netCDF for one model/scenario, written by grow_seas.ipynb."""
    return f"{cache_dir}/phen_{scenario}_{model}.nc"


def load_growing_season(models, scenario="ctl", cache_dir=GS_CACHE_DIR):
    """Load the cached per-model phenology Datasets (onset=gst_start_month,
    offset=gst_end_month, gst_type) written by grow_seas.ipynb. Raises a clear
    FileNotFoundError (rather than xarray's less obvious one) if grow_seas.ipynb
    has not been run yet."""
    import os
    out = []
    for mm in models:
        p = gs_cache_path(scenario, mm, cache_dir)
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Growing-season cache not found: {p}\n"
                f"Run grow_seas.ipynb first (it writes this file for every model/scenario)."
            )
        out.append(xr.open_dataset(p))
    return out
