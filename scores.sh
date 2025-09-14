#! /bin/bash

for a in output/*.msf
do 
	/home/ubuntu/mnt/BAliBASE_R9/src/bali_score ~/mnt/BAliBASE/RV100/BBA0001.xml $a > eca
	echo $a | cut -d. -f2-4
	cat eca | grep SP | cut -d' ' -f3
	cat eca | grep TC | cut -d' ' -f3

	/home/ubuntu/mnt/bali_score_src/bali_score ~/mnt/BAliBASE/RV100/BBA0001.xml $a > eca
	cat eca | grep CS | cut -d' ' -f3
done
