#!/bin/bash

####################################################################################
# Script to make all points mask	
####################################################################################

lonmax=179
latmin=-89

gmt xyz2grd all_points.xyz -An -I1 -Rd -fg -Gjunk.grd
gmt grd2xyz junk.grd | awk '{if ($3!="NaN") {print $1, $2, 1}else{print $1, $2, 0}}' | sort -k2,2rn -k1,1n | \
	awk '{if ($2>='$latmin' && $1<='$lonmax'){print $3}}' > all_points_mask_new.z 

rm junk.grd