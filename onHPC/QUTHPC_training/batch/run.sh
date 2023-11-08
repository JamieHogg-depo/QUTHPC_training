#!/bin/bash

base_folder='QUTHPC_training/batch'

# Model
for model_spec in asra1 asra2 ## Name of file that runs the model
do
	
	# Sex
	for sex in Persons Female Male
	do
	
		# Condition
		for condition in Asthma CHD
		do
	
			# create the unique file name
			Rfile=$model_spec'_'$condition'_'$sex
	
			# get the current date
			cur_date=$(date +%Y%m%d%H%M)
		
			# create directories
			mkdir -p $base_folder/sub_src/$cur_date
			mkdir -p $base_folder/outputs/lyra/$cur_date/errors
			mkdir -p $base_folder/outputs/lyra/$cur_date/out
			mkdir -p $base_folder/outputs/$condition/$sex
			
			# Set loop output file with full directory
			loop_output_file=$base_folder'/outputs/'$condition'/'$sex'/'$Rfile

			# create the unique .sub script files
			file=$Rfile'.sub'

			# paste the commands in the .sub scripts
			cat > $base_folder/sub_src/$cur_date/$file <<EOF
#!/bin/bash -l
#PBS -N $Rfile
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=0:20:00
#PBS -e $base_folder/outputs/lyra/$cur_date/errors/$Rfile
#PBS -o $base_folder/outputs/lyra/$cur_date/out/$Rfile

module load r/4.0.3-foss-2020b
module load gdal/3.2.1-foss-2020b

R -e ".libPaths('r_lib');
base_folder='$base_folder';
loop_output_file='$loop_output_file';
model_spec='$model_spec';
sex='$sex';
condition='$condition';
cur_date='$cur_date';
niter=800;
nburnin=400;
thin=1;
nchains=4;
source('$base_folder/ms.R');"
EOF

			# run each script
			qsub $base_folder/sub_src/$cur_date/$file

		done
	done
done
