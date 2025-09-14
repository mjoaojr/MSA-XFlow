#/bin/bash

FRMW_DIR=~/mnt/framework

frmw_full=~/mnt/framework/full/src/frmw_full
frmw_quick=~/mnt/framework/quick/src/frmw_quick
frmw_kmers=~/mnt/framework/kmers/src/frmw_kmers
frmw_lcs=~/mnt/framework/lcs/src/frmw_lcs
frmw_upgma=~/mnt/framework/upgma/src/frmw_upgma
frmw_slmin=~/mnt/framework/slmin/src/frmw_slmin
frmw_slmax=~/mnt/framework/slmax/src/frmw_slmax
frmw_nj=~/mnt/framework/nj/src/frmw_nj
frmw_align=~/mnt/framework/align/src/frmw_align
mclustalw=~/mnt/mclustalw/src/clustalw2

cd work 
if [ -e saida ]
then
	rm saida
fi

for eca in $1/*.in_tfa 
do
	arq=$eca
#	echo $arq
	base=$(basename $arq .in_tfa)

	for score in QUICKTREE LCS FULL KMERS
#	for score in QUICKTREE LCS FULL KMERS PROBA
	do
		echo "****************** $score *****************"
		if [ "$score" = 'FULL' ]
		then
			$frmw_full -ALIGN -INFILE=$arq >> saida
			ss='full'
		else
			if [ $score = 'QUICKTREE' ]
			then
				$frmw_quick -ALIGN -INFILE=$arq >> saida
				ss='quick'
			else
				if [ $score = 'KMERS' ]
				then
					$frmw_kmers -ALIGN -INFILE=$arq >> saida
					ss='kmers'
				else
					$frmw_lcs -ALIGN -INFILE=$arq >> saida
					ss='lcs'
				fi
			fi
		fi

#		echo $base.$ss.ss
		mv $1/$base.$ss.ss .

		for cluster in NJ UPGMA SLMIN SLMAX
		do
			echo "================== $score $cluster ================="
			if [ $cluster = 'NJ' ]
			then
				$frmw_nj -ALIGN -INFILE=$arq -SSFILE=$base.$ss.ss >> saida
				gt='nj'
			else
				if [ $cluster = 'UPGMA' ]
				then
					$frmw_upgma -ALIGN -INFILE=$arq -SSFILE=$base.$ss.ss >> saida
					gt='upgma'
				else
					if [ $cluster = 'SLMIN' ]
					then
						$frmw_slmin -ALIGN -INFILE=$arq -SSFILE=$base.$ss.ss >> saida
						gt='slmin'
					else
						$frmw_slmax -ALIGN -INFILE=$arq -SSFILE=$base.$ss.ss >> saida
						gt='slmax'
					fi
				fi
			fi
			
#			echo $base.$ss.$gt.dnd
			mv $1/$base.dnd ./$base.$ss.$gt.dnd
 			
			$frmw_align -ALIGN -CLUSTERING=$cluster -INFILE=$arq -SSFILE=$base.$ss.ss -USETREE=$base.$ss.$gt.dnd -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$base.$ss.$gt.msf -OUTORDER=INPUT > saida
#			echo $base.$ss.$gt.msf

#			if [ $score = 'FULL' ]
#			then
#				$mclustalw -ALIGN -CLUSTERING=$cluster -INFILE=$arq -CLUSTERING=$cluster -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$base.$ss.$gt-mclustal.msf -OUTORDER=INPUT > saida
#				mv $1/$base.dnd ./$base.full.$gt-mclustal.dnd
#			else
#				$mclustalw -ALIGN -$score -CLUSTERING=$cluster -INFILE=$arq -MAXDIV=0 -OUTPUT=GCG -OUTFILE=$base.$ss.$gt-mclustal.msf -OUTORDER=INPUT > saida
#				mv $1/$base.dnd ./$base.$ss.$gt-mclustal.dnd
#			fi
#
#			diff $base.$ss.$gt.msf $base.$ss.$gt-mclustal.msf
#			diff $base.$ss.$gt.dnd $base.$ss.$gt-mclustal.dnd

	#               /usr/bin/time -p ./testa.sh $1 -$score $cluster
	#		./testa.sh $1 -$score $cluster
	#               cat spscore
	#		cat myspscore
	#               /usr/bin/time -p ./testa.sh $1 -$score $cluster -CONSISTENCY
	#               cat spscore
	#               /usr/bin/time -p ./testa.sh $1 -$score $cluster -ITERATION=muscle -NUMITER=100
	#               cat spscore
	#               /usr/bin/time -p ./testa.sh $1 -$score $cluster -CONSISTENCY -ITERATION=muscle -NUMITER=100
	#               cat spscore
		done
	done
done
