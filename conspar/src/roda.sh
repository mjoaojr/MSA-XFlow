#!/bin/bash

for i in `seq 1` 
do
	/usr/bin/time -p ./frmw_conspar-$1 -ALIGN -INFILE=pf00005-split_$2-$3-$4.in_tfa -NTHR=$5 > eca
	grep Tempo eca
done
