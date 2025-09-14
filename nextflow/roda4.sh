#!/bin/bash

if [ ! -f ../input/enable/pf00005-split_$2-$3-$4.in_tfa ] 
then
	rm ../input/enable/*
	cp ../input/seqs/pf00005-split_$2-$3-$4.in_tfa ../input/enable/
fi

if [ -f tempo ] 
then
	rm tempo
fi

for i in `seq $5 -1 0`
do
	j=$(( $1 - $i ))
	echo ============================================= $j $i ================================================ >> tempo
	rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-conspar.nf -c 4cpu.config -with-timeline time$1conspar_$2-$3-$4.$j-$i.html --in ~/local/framework/input/enable/ --n $i 2>> tempo
	rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-conspar.nf -c 4cpu.config -with-timeline time$1conspar_$2-$3-$4.$j-$i.html --in ~/local/framework/input/enable/ --n $i 2>> tempo
	rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-conspar.nf -c 4cpu.config -with-timeline time$1conspar_$2-$3-$4.$j-$i.html --in ~/local/framework/input/enable/ --n $i 2>> tempo
done
