#!/bin/bash


############################
# Resample at the same interval as PREM
############################

awk '{print $1}' ../PREM.dat > knotfile
scalingfactor010="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta=0.10.txt"
scalingfactor009="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta=0.09.txt"
scalingfactor008="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta=0.08.txt"
scalingfactor007="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta=0.07.txt"
scalingfactor006="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta=0.06.txt"

scalingfactor010_2="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta_2=0.10.txt"
scalingfactor009_2="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta_2=0.09.txt"
scalingfactor008_2="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta_2=0.08.txt"
scalingfactor007_2="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta_2=0.07.txt"
scalingfactor006_2="/space/mch74/gmt/RESIDUAL_DEPTH_SH/02_figures/SH_power_spectra_analytical/output/scalingfactor_Beta_2=0.06.txt"

echo ""
pwd
rm -f FINAL_density.dat
touch FINAL_density.dat

for model in J362D28 S20A-ISO S20RTS SAW24B16 SB4_L18 TX2002; do
	echo "Density Model: $model"
	awk '{print $1, $2}' density_${model}.dat | gmt sample1d -Fa -Ar -Tknotfile > density_${model}_resampled.dat
	awk '{print $1, $2}' swave_${model}.dat | gmt sample1d -Fa -Ar -Tknotfile > swave_${model}_resampled.dat
	paste density_${model}_resampled.dat swave_${model}_resampled.dat | awk '{print $1, $2, $4, $2/$4,(($2/100)*3330)}' > FINAL_${model}.dat
	awk '{print $1, $2}' $scalingfactor010 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_010_resampled.dat
	awk '{print $1, $2}' $scalingfactor009 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_009_resampled.dat
	awk '{print $1, $2}' $scalingfactor008 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_008_resampled.dat
	awk '{print $1, $2}' $scalingfactor007 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_007_resampled.dat
	awk '{print $1, $2}' $scalingfactor006 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_006_resampled.dat

	awk '{print $1, $2}' $scalingfactor010_2 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_2_010_resampled.dat
	awk '{print $1, $2}' $scalingfactor009_2 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_2_009_resampled.dat
	awk '{print $1, $2}' $scalingfactor008_2 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_2_008_resampled.dat
	awk '{print $1, $2}' $scalingfactor007_2 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_2_007_resampled.dat
	awk '{print $1, $2}' $scalingfactor006_2 | gmt sample1d -Fa -Ar -Tknotfile > scalingfactor_2_006_resampled.dat
done

############################
# collate all tomographic_models into single file, one for delta rho / rho, one for delta Vs/ Vs, one for the ratio
############################

paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat | awk '{print $1, $2, $7, $12, $17, $22, $27}' > FINAL_density.dat
paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat | awk '{print $1, $3, $8, $13, $18, $23, $28}' > FINAL_swave.dat
paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat | awk '{print $1, $4, $9, $14, $19, $24, $29}' > FINAL_ratio.dat
paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat | awk '{print $1, $5, $10, $15, $20, $25, $30}' > FINAL_density2.dat

# OPTION 5 - using the mean of the scaling factors
var="10"
paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat FINAL_ratio_average.dat | awk -v var=$var '{print $1, $3*var*3.33*$32, $8*var*3.33*$32, $13*var*3.33*$32, $18*var*3.33*$32, $23*var*3.33*$32, $28*var*3.33*$32}' > FINAL_density_variablemean.dat

# OPTION 6,7,8,9,10 - using a cosine approximation. Try different beta factors
var="10"
for Beta in 010 009 008 007 006; do 
	paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat scalingfactor_${Beta}_resampled.dat | awk -v var=$var '{print $1, $3*var*3.33*$32, $8*var*3.33*$32, $13*var*3.33*$32, $18*var*3.33*$32, $23*var*3.33*$32, $28*var*3.33*$32}' > FINAL_density_variablemean_${Beta}_cosine.dat
done

for Beta in 010 009 008 007 006; do 
	paste FINAL_J362D28.dat FINAL_S20A-ISO.dat FINAL_S20RTS.dat FINAL_SAW24B16.dat FINAL_SB4_L18.dat FINAL_TX2002.dat scalingfactor_2_${Beta}_resampled.dat | awk -v var=$var '{print $1, $3*var*3.33*$32, $8*var*3.33*$32, $13*var*3.33*$32, $18*var*3.33*$32, $23*var*3.33*$32, $28*var*3.33*$32}' > FINAL_density_variablemean_2_${Beta}_cosine.dat
done

################################################
# find mean, min, max
################################################

echo "0 12" > top.temp
echo "2900 12" > tail.temp

###################
# find average
###################
# density
awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_density.dat > FINAL_density_average.dat

# swave
awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_swave.dat > FINAL_swave_average.dat
cat top.temp FINAL_swave_average.dat tail.temp > ../FINAL_swave_average.dat
# ratio
awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_ratio.dat > FINAL_ratio_average.dat


# density 2
awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_density2.dat > FINAL_density2_average.dat
cat top.temp FINAL_density2_average.dat tail.temp > ../FINAL_density2_average.dat

# density - variable mean
awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_density_variablemean.dat > FINAL_density_variablemean_average.dat
cat top.temp FINAL_density_variablemean_average.dat tail.temp > ../FINAL_density_variablemean_average.dat

# density - cosine

for Beta in 010 009 008 007 006; do 
	# generate average
	awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_density_variablemean_${Beta}_cosine.dat > FINAL_density_variablemean_${Beta}_cosine_average.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_${Beta}_cosine_average.dat tail.temp > ../FINAL_density_variablemean_${Beta}_cosine_average.dat
done

for Beta in 010 009 008 007 006; do 
	# generate average
	awk '{print $1, ($2+$3+$4+$5+$6+$7)/6}' FINAL_density_variablemean_2_${Beta}_cosine.dat > FINAL_density_variablemean_2_${Beta}_cosine_average.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_2_${Beta}_cosine_average.dat tail.temp > ../FINAL_density_variablemean_2_${Beta}_cosine_average.dat
done

###################
# find min
###################
# density
awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_density.dat > FINAL_density_min.dat

# swave
awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_swave.dat > FINAL_swave_min.dat
cat top.temp FINAL_swave_min.dat tail.temp > ../FINAL_swave_min.dat
# ratio
awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_ratio.dat > FINAL_ratio_min.dat

# density 2
awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_density2.dat > FINAL_density2_min.dat
cat top.temp FINAL_density2_min.dat tail.temp > ../FINAL_density2_min.dat

# density - variable mean
awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_density_variablemean.dat > FINAL_density_variablemean_min.dat
cat top.temp FINAL_density_variablemean_min.dat tail.temp > ../FINAL_density_variablemean_min.dat

# density - cosine

for Beta in 010 009 008 007 006; do 
	# generate min
	awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_density_variablemean_${Beta}_cosine.dat > FINAL_density_variablemean_${Beta}_cosine_min.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_${Beta}_cosine_min.dat tail.temp > ../FINAL_density_variablemean_${Beta}_cosine_min.dat
done

for Beta in 010 009 008 007 006; do 
	# generate min
	awk '{max=$2;for(i=2;i<=NF;i++){if($i > max) max = $i}print $1, max}' FINAL_density_variablemean_2_${Beta}_cosine.dat > FINAL_density_variablemean_2_${Beta}_cosine_min.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_2_${Beta}_cosine_min.dat tail.temp > ../FINAL_density_variablemean_2_${Beta}_cosine_min.dat
done

###################
# find max
###################
# density
awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_density.dat > FINAL_density_max.dat

# swave
awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_swave.dat > FINAL_swave_max.dat
cat top.temp FINAL_swave_max.dat tail.temp > ../FINAL_swave_max.dat
# ratio
awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_ratio.dat > FINAL_ratio_max.dat

# density 2
awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_density2.dat > FINAL_density2_max.dat
cat top.temp FINAL_density2_max.dat tail.temp > ../FINAL_density2_max.dat

# density - variable mean
awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_density_variablemean.dat > FINAL_density_variablemean_max.dat
cat top.temp FINAL_density_variablemean_max.dat tail.temp > ../FINAL_density_variablemean_max.dat

# density - cosine

for Beta in 010 009 008 007 006; do 
	# generate max
	awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_density_variablemean_${Beta}_cosine.dat > FINAL_density_variablemean_${Beta}_cosine_max.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_${Beta}_cosine_max.dat tail.temp > ../FINAL_density_variablemean_${Beta}_cosine_max.dat
done

for Beta in 010 009 008 007 006; do 
	# generate max
	awk '{min=$2;for(i=2;i<=NF;i++){if($i < min) min = $i}print $1, min}' FINAL_density_variablemean_2_${Beta}_cosine.dat > FINAL_density_variablemean_2_${Beta}_cosine_max.dat
	# generate final curve
	cat top.temp FINAL_density_variablemean_2_${Beta}_cosine_max.dat tail.temp > ../FINAL_density_variablemean_2_${Beta}_cosine_max.dat
done

rm top.temp
rm tail.temp

############################
# top and tail and move folder
############################
# to do



cat top.temp tail.temp


echo ""
echo "Finished collating data"
echo ""

############################
# Plot tomographic_models
############################

echo ""
echo "Begin Plotting"
echo ""

# -------------------------------------
# --------- SET GMT DEFAULTS ----------
# -------------------------------------

gmt gmtset MAP_FRAME_TYPE = plain
gmt gmtset PS_MEDIA = a0
gmt gmtset FONT_ANNOT_PRIMARY = 8
gmt gmtset MAP_ANNOT_OFFSET = 0.1
gmt gmtset FONT_LABEL = 12
gmt gmtset MAP_TICK_PEN_PRIMARY = 0.5
gmt gmtset MAP_FRAME_PEN = 1
gmt gmtset MAP_TICK_PEN = 1

psfile="temp.temp"

############################
# DELTA DENISTY / DENSITY
############################

rgn="-R0/3000/0.0/0.4"
proj="-JX10/6"

echo "Plot - delta density / density"

gmt psbasemap $rgn $proj -Ba0 -K -Y60c > $psfile
# plot all profiles

for model in J362D28 S20A-ISO S20RTS SAW24B16 SB4_L18 TX2002; do
	awk '{print $1, $2}' FINAL_${model}.dat | gmt psxy $rgn $proj -W0.5,black -O -K >> $psfile
done
#awk '{print $1, $2}' FINAL_density_average.dat | gmt psxy $rgn $proj -W1,blue -O -K >> $psfile
#awk '{print $1, $2}' FINAL_density_min.dat | gmt psxy $rgn $proj -W1,blue -O -K >> $psfile
#awk '{print $1, $2}' FINAL_density_max.dat | gmt psxy $rgn $proj -W1,blue -O -K >> $psfile

echo "a" | gmt pstext $rgn $proj  -D12p/-12p -F+f22p,Helvetica,black+jTL+cTL -Gwhite -TO -W0.5,black -N -O -K >> $psfile

gmt psbasemap $rgn $proj -Ba500f100:"Depth, km":/a0.1f0.05:"@~\144@~@~\162@~ / @~\162@~ (%)":WeSn -O -K >> $psfile

############################
# DELTA DENISTY 
############################

rgn="-R0/3000/0.0/12"
proj="-JX10/6"

echo "Plot - density"

gmt psbasemap $rgn $proj -Y-7 -Ba0 -K -O >> $psfile

# plot all profiles

for model in J362D28 S20A-ISO S20RTS SAW24B16 SB4_L18 TX2002; do
	awk '{print $1, $2}' density2_${model}.dat | gmt psxy $rgn $proj -W0.5,gray -O -K >> $psfile
done

awk '{print $1, $2}' FINAL_density2_average.dat | gmt psxy $rgn $proj -W2,purple -O -K >> $psfile
awk '{print $1, $2}' FINAL_density2_min.dat | gmt psxy $rgn $proj -W2,purple -O -K >> $psfile
awk '{print $1, $2}' FINAL_density2_max.dat | gmt psxy $rgn $proj -W2,purple -O -K >> $psfile


# this is what we want to replicate
awk '{print $2, $1}' ../average_density_anomalies.dat | gmt psxy $rgn $proj -W0.5,black -O -K >> $psfile
awk '{print $2, $1}' ../maximum_density_anomalies.dat | gmt psxy $rgn $proj -W0.5,black,-- -O -K >> $psfile
awk '{print $2, $1}' ../minimum_density_anomalies.dat | gmt psxy $rgn $proj -W0.5,black,-- -O -K >> $psfile

#awk '{print $1, $2}' FINAL_density_average.dat | gmt psxy $rgn $proj -W1,navy -O -K >> $psfile
echo "b" | gmt pstext $rgn $proj  -D12p/-12p -F+f22p,Helvetica,black+jTL+cTL -Gwhite -TO -W0.5,black -N -O -K >> $psfile


gmt psbasemap $rgn $proj -Ba500f100:"Depth, km":/a4f1:"@~\144@~@~\162@~ ":WeSn -O -K >> $psfile
############################
# SHEAR WAVE
############################

rgn="-R0/3000/0.0/3"
proj="-JX10/6"

echo "Plot - shear wave"

gmt psbasemap $rgn $proj -Y-7 -Ba0 -K -O >> $psfile
for model in J362D28 S20A-ISO S20RTS SAW24B16 SB4_L18 TX2002; do
	awk '{print $1, $3}' FINAL_${model}.dat | gmt psxy $rgn $proj -W0.5,gray -O -K >> $psfile
done
awk '{print $1, $2}' FINAL_swave_average.dat | gmt psxy $rgn $proj -W1,red -O -K >> $psfile
awk '{print $1, $2}' FINAL_swave_min.dat | gmt psxy $rgn $proj -W1,red -O -K >> $psfile
awk '{print $1, $2}' FINAL_swave_max.dat | gmt psxy $rgn $proj -W1,red -O -K >> $psfile

gmt psbasemap $rgn $proj -Ba500f100:"Depth, km":/a0.5f0.1:"@~\144@~V@-S@- / V@-S@-":WeSn -O -K >> $psfile
echo "c" | gmt pstext $rgn $proj  -D12p/-12p -F+f22p,Helvetica,black+jTL+cTL -Gwhite -TO -W0.5,black -N -O -K >> $psfile

############################
# Velocity–density scaling factor
############################

rgn="-R0/3000/0/0.3"
proj="-JX10/6"

echo "Plot - velocity-density scaling factor"

gmt psbasemap $rgn $proj -Y-7 -Ba0 -K -O >> $psfile
for model in J362D28 S20A-ISO S20RTS SAW24B16 SB4_L18 TX2002; do
	awk '{print $1, $4}' FINAL_${model}.dat | gmt psxy $rgn $proj -W0.5,gray -O -K >> $psfile
done
awk '{print $1, $2}' FINAL_ratio_average.dat | gmt psxy $rgn $proj -W1,hotpink -O -K >> $psfile
echo "d" | gmt pstext $rgn $proj  -D12p/-12p -F+f22p,Helvetica,black+jTL+cTL -Gwhite -TO -W0.5,black -N -O -K >> $psfile
gmt psbasemap $rgn $proj -Ba500f100:"Depth, km":/a0.1f0.05:"velocity-density scaling factor":WeSn -O >> $psfile

echo "Plot - convert to jpg"
jpgfile=$(echo $psfile | sed s/"\.ps"/"\.jpg"/)
gmt psconvert -Tj -A+m0.05c -P -E300 -Z $psfile

rm -f gmt.conf
rm -f gmt.history
