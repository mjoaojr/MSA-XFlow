#!/bin/bash

if [ ! -f ../input/enable/pf00005-split_$1-$2-$3.in_tfa ] 
then
	rm ../input/enable/*
	cp ../input/seqs/pf00005-split_$1-$2-$3.in_tfa ../input/enable/
fi

if [ -f tempo ] 
then
	rm tempo
fi

echo ============================================= $1 $2 $3 ================================================ >> tempo
rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-cons-regr.nf -c 1cpu.config -with-timeline timecons-regr_$1-$2-$3.html --in ~/local/framework/input/enable/ 2>> tempo
rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-cons-regr.nf -c 1cpu.config -with-timeline timecons-regr_$1-$2-$3.html --in ~/local/framework/input/enable/ 2>> tempo
rm -rf work/*; /usr/bin/time -p nextflow -q run frmw-cons-regr.nf -c 1cpu.config -with-timeline timecons-regr_$1-$2-$3.html --in ~/local/framework/input/enable/ 2>> tempo
