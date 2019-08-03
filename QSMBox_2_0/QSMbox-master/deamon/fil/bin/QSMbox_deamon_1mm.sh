#!/bin/bash

#ptbsfile=ptbs_cme_002_def_msdi2
ptbsfile=ptbs_ppd_001_def_msdi_1mm
forcerun=0

d0=$PWD
TS=$(date +"%Y%m%d_%H%M%S")

msdi2deamondir=$(dirname $(jpath $0))
cd $msdi2deamondir/..
d1=$PWD

cd $d0
datadirs=$(jdirs)

if [ -n "$datadirs" ]; then
	for i in $datadirs; do echo ++${i}
		if [ -f "$d0/$i/${ptbsfile}.m" ]; then
			echo PTBS file found in $i, archive it and import preset file from ${d1}.
			mv -v "$d0/$i/${ptbsfile}.m" "$d0/$i/ARCH_${ptbsfile}_${TS}.m"
			cp -v $d1/${ptbsfile}.m $d0/$i/
			prev_flag=1
		else
			echo PTBS file not found in $i, import from ${d1}.
			cp -v $d1/${ptbsfile}.m $d0/$i/
			prev_flag=0
		fi
		
		cd $d0/$i

		if ([ "$prev_flag" == "1" ] && [ ! -f "$d0/$i/msdiout.txt" ]) || [ "$prev_flag" == "0" ] || [ "$forcerun" == "1" ]; then
			jmatcml $ptbsfile
		else
			echo "A previous PTBS file with the same name was found in the working dir, $d0/$i/msdiout.txt exists and forcerun set to 0 => Do nothing."
		fi

		echo
	done
else
	echo No data dirs found in ${d0}.
fi

