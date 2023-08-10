#!/bin/bash

#for var in 'nME' 'wME' 'uwME'
for var in 'uwME'
#for var in 'NewWeights'
do

# cur_data
Rfile='YLDMeasurementError'
#Rfile='Smoking_MrP'
Rfile_c=$Rfile'_'$var
cur_date=$(date +%Y%m%d)
	
# create directories
mkdir -p WADOH/sub_src/$cur_date
mkdir -p WADOH/outputs/$cur_date/lyra_errors
mkdir -p WADOH/outputs/$cur_date/lyra_out
mkdir -p WADOH/outputs/$cur_date/r

# create the unique .sub script files
file=$Rfile_c'.sub'

# paste the commands in the .sub scripts
cat>WADOH/sub_src/$cur_date/$file<<EOF
#!/bin/bash -l
#PBS -N $Rfile_c
#PBS -l ncpus=1
#PBS -l mem=60GB
#PBS -l walltime=50:00:00
#PBS -e WADOH/outputs/$cur_date/lyra_errors/$Rfile_c
#PBS -o WADOH/outputs/$cur_date/lyra_out/$Rfile_c

module load r/4.0.3-foss-2020b
module load gdal/3.2.1-foss-2020b

R -e ".libPaths('r_lib');
var='$var'
Rfile='$Rfile_c';
cur_date='$cur_date';
niter=80000;
nburnin=40000;
thin=20;
nchains=4;
source('WADOH/master.R');"
EOF

	# run each script
		qsub WADOH/sub_src/$cur_date/$file

done
			