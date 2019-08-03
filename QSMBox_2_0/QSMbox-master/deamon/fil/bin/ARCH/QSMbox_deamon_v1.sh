#!/bin/bash

msdi2deamondir=$(dirname $(jpath $0))
ptbsfile=ptbs_cme_002_def_msdi2

#recondir=$1
recondir=$PWD

forcerun=0 #$2

cd $recondir
datadirs=$(jdirs)

if [ -n "$datadirs" ]; then
	for i in $datadirs; do echo ++${i}
		if [ -f "$recondir/$i/${ptbsfile}.m" ]; then
			if [ "$forcerun" == "1" ]; then
				echo PTBS file found in $i, archive it and import new file from $msdi2deamondir.
				mv -v "$recondir/$i/${ptbsfile}.m" "$recondir/$i/bku_${ptbsfile}.m"
				cp -v $msdi2deamondir/../${ptbsfile}.m $recondir/$i/
			else
				echo PTBS file found in $i, do not import.
			fi
		else
			echo PTBS file not found in $i, import from $msdi2deamondir.
			cp -v $msdi2deamondir/../${ptbsfile}.m $recondir/$i/
		fi
		
		cd $recondir/$i

		if [ "$forcerun" == "1" ]; then
			qsmfiles=$(ls qsm_INTEGRAL_*.nii* 2>/dev/null)
			if [ -n "$qsmfiles" ]; then
				echo Archive QSM_INTEGRAL files.
				for j in $qsmfiles; do echo $j
					mv -v $j bku_${j}
				done
			fi
		fi

		if [ ! -f "$recondir/$i/msdiout.txt" ] || [ "$forcerun" == "1" ]; then
			jmatcml $ptbsfile
		else
			echo "$recondir/$i/msdiout.txt exists and \$forcerun=0, do nothing."
		fi
	done
else
	echo No data dirs found in ${recondir}.
fi
