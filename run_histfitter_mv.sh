#!/bin/bash

mv_list="0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 
0.105 0.11 0.115 0.12 0.125 0.13 0.135 0.14 0.145 0.15 0.155 0.16 0.165 0.17 0.175 0.18 0.185 0.19 0.195 0.2 0.205 0.21 
0.215 0.22 0.225 0.23 0.235 0.24 0.245 0.25 0.255 0.26 0.265 0.27 0.275 0.28 0.285 0.29 0.295 0.3 0.305 0.31 0.315 0.32 
0.325 0.33 0.335 0.34 0.345 0.35 0.355 0.36 0.365 0.37 0.375 0.38 0.385 0.39 0.395 0.4 0.405 0.41 0.415 0.42 0.425 0.43 
0.435 0.44 0.445 0.45 0.455 0.46 0.465 0.47 0.475 0.48 0.485 0.49 0.495 0.5 0.505 0.51 0.515 0.52 0.525 0.53 0.535 0.54 0.545"

mv_list0="0.005 0.01 0.015"

rm 000.dat 111.dat
for mv in $mv_list
do
 logfile=./results/ldm$mv.log
 echo $mv $logfile
 HistFitter.py  -t -w -D "before" -f -F excl -l  -V ./FGDConfig_v1.py -L DEBUG -u $mv > $logfile
 explimit=`grep median $logfile | awk '{print $6}'`
 datalimit=`grep "computed upper limit" $logfile | awk '{print $8}'`
 echo $mv  $explimit  $datalimit
 echo $mv  $explimit >> 000.dat
 echo $mv  $datalimit >> 111.dat
done

echo Done...
