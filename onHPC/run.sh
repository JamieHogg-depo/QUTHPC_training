#!/bin/bash

base_folder='QUTHPC_training/onHPC'

# Model
for model_spec in model1 #model2
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
	cur_date=$(date +%Y%m%d)
		
	# create directories for current date
	mkdir -p $base_folder/sub_src/$cur_date
	mkdir -p $base_folder/outputs/$cur_date/lyra_errors
	mkdir -p $base_folder/outputs/$cur_date/lyra_out
	mkdir -p $base_folder/outputs/$cur_date/r

	# create the unique .sub script files
	file=$Rfile'.sub'

	# paste the commands in the .sub scripts
	cat>$base_folder/sub_src/$cur_date/$file<<EOF
#!/bin/bash -l
#PBS -N $Rfile
#PBS -l ncpus=1
#PBS -l mem=2GB
#PBS -l walltime=1:00:00
#PBS -e $base_folder/outputs/$cur_date/lyra_errors/$Rfile
#PBS -o $base_folder/outputs/$cur_date/lyra_out/$Rfile

module load r/4.0.3-foss-2020b
module load gdal/3.2.1-foss-2020b

R -e ".libPaths('r_lib');
base_folder='$base_folder';
model_spec='$model_spec';
sex='$sex';
condition='$condition';
Rfile='$Rfile';
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
			