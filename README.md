This repository contains the codes for the replication of the results presented in the following manuscript under review in AGU Earth's Future Journal:
"Ecohydrological regime modulates the thermal and water-cycle response to large-scale afforestation in a high CO2 climate", by Mastropierro et al., (2026).

The repository is organized as follows.

```bash
/Afforestation_LUMIP_biophys
│   README.md
│   environment.yml
└───Code
        lumip_utils.py
        grow_seas.ipynb
        mlr_coomp.ipynb
        mlr_proc.ipynb
        rough_tran_scale.ipynb
        seb_grow_seas.ipynb
        water_avail_final.ipynb
```

**environment.yml** is the conda environment containing the libraries used to perform the analysis.

The folder **Code** contains the .ipynb files, used to process the data and generate the figures of the manuscript, and the lumip_utils.py with some common utility functions used by the notebooks. 
