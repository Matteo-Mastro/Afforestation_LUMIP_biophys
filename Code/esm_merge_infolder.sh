#!/bin/bash
sudo mount -t drvfs F: /mnt/f

cd "/mnt/f/Data/LUMIP"
base="/mnt/f/Data/LUMIP"

scenario=(ssp370 ssp370-ssp126Lu)
model=(UKESM1-0-LL MPI-ESM1-2-LR)

vars=(evspsblveg_Lmon tran_Lmon) # va_Amon ua_Amon

for mm in ${model[@]}; do  		## Loop over all models
		
		for ss in ${scenario[@]}; do		## Loop over all scenarios
		# realiz=$(ls ${base}/${ss}/${mm})		## Take the realizations' folders inside each model

			# for rr in ${realiz[@]}; do			## Loop over all the realizations (overcome different forcings prob)
				in_dir="${base}/${mm}"

				for vv in ${vars[@]}; do
				f=$(ls ${in_dir}/${vv}_*)
				
				if [ ${ss} = historical ]; then

					# if exist multiple files in time domain (with the same name)
					if [ `ls -1 $f | wc -l ` -gt 1 ]; then			
					echo "${mm} has multiple files in ${ss}, merging ${vv}..."				
					cdo mergetime ${f[@]} ${in_dir}/${vv}_${mm}_${ss}_${rr}_gn_185001-201412.nc
            		else								# is files are already merged in time
					echo "${mm} already merged in ${ss}, skipping..."
					fi

				else

					# if exist multiple files in time domain (with the same name)
					if [ `ls -1 $f | wc -l ` -gt 1 ]; then			
					echo "${mm} has multiple files in ${ss}, merging ${vv}..."	
					cdo mergetime ${f[@]} ${in_dir}/${vv}_${mm}_${ss}_gn_201501-210012.nc
					else								# is files are already merged in time
					echo "${mm} already merged in ${ss}, skipping..."
					fi
				fi 

				done # vv

			# done # rr

		done # ss

done # ss

exit


