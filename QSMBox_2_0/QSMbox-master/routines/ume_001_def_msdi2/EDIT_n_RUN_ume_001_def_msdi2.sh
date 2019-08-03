# Instructions
# ============
#
# 1. All sM*.nii and sP*.nii files must be copied/moved to this directory.
#
#  Nb. This assumes you converted nifti files from dicom using SPM's
#  'Dicom Import' utility
#
#
# 2. Import dicom file (this will be used to read in a few scan params).
#
#  mkdir dicom
#
#  cp sample dicom file to the 'dicom' dir
#
#  If all your dicom files are in one directory, you can search for it
#  printing the dicom 'Series Description' field for all files as follows:
#
#  for i in *dcm; do echo $i; mri_probedicom --i ${i} --t 0008 103e; done
#
#  Nb. mri_probedicom comes with FreeSurfer
#
#
# 3. Edit input arguments (below), and run this file:
#
#
# Created by Julio Acosta-Cabronero


################################## INPUT ######################################
TE1=4.92e-3
TE2=9.84e-3 
TE3=19.70e-3               # Add more TEs if necessary
deltaTE=4.92e-3            # deltaTE=TE2-TE1
TEbrainmask=$TE3           # TE~20ms (3T) is usually a good choice for step #1
TEall="$TE1 $TE2 $TE3"     # Add more TEs if necessary
number_of_coils=32         # # of receive channels
phase_scaling_factor=4096  # Phase ranging (-4096,4096) is Siemens convention
###############################################################################


###############################################################################

echo
echo 1/14. Prep for brain mask calculation

d0=$PWD

rm -rf $d0/brainmask

mkdir $d0/brainmask

rm -f ptb_*.txt 

echo $TEbrainmask > ptb_TE.txt

echo $phase_scaling_factor > ptb_scaling_factor.txt

cp -R dicom brainmask/
cp ptb_TE.txt brainmask/
cp ptb_scaling_factor.txt brainmask/
cp ptbs_ume_001_step1.m brainmask/
cp jmatcml brainmask/

###############################################################################

echo
echo 2/14. Copy 3rd echo data to brainmask directory

for i in `seq 1 $number_of_coils`; do echo $i
	for j in M P; do 
		echo $j
		cp s${j}${i}_3_* brainmask/
	done
done

###############################################################################

echo
echo "3/14. Merge magn_orig (3rd echo uncombined)"

cd $d0/brainmask

fslmerge -t magn_orig sM*

###############################################################################

echo
echo "4/14. Merge phase_orig (3rd echo uncombined)"

fslmerge -t phase_orig sP*

rm -f sM* sP*

###############################################################################

echo
echo "5/14. Calculate brainmask (step #1)"

./jmatcml ptbs_ume_001_step1

###############################################################################

echo
echo 6/14. Move all nifit files to channel-specific directories

for i in `seq 1 $number_of_coils`; do 
	for j in M P; do 
		N=$(./jnums $i $i 2)
		echo + $N
		if [ "$j" == "M" ]; then
			mkdir ch${N}
		fi
		mv s${j}${i}_* ch${N}/
	done
done

###############################################################################

echo
echo "7/14. Merge multiple echoes into 4D files (one per coil element)"

for i in ch*; do echo $i
	cd $d0/$i
	fslmerge -t magn_orig sM*
	fslmerge -t phase_orig sP*
done	

###############################################################################

echo
echo "8/14. Archive original nifti files in orig (directory)"

cd $d0

mkdir orig

mv ch*/*nii orig/

###############################################################################

echo
echo "9/14. Echo combination (step #2)"

rm -f ptb_TE.txt

echo $TEall > ptb_TE.txt

for i in ch*; do echo $i

	cd $d0
	cp -R dicom $i/
	cp ptb_TE.txt $i/
	cp ptb_scaling_factor.txt $i/
	cp ptbs_ume_001_step2.m $i/
	cp brainmask/roi.nii $i/
	cp jmatcml $i/

	cd $d0/$i
	./jmatcml ptbs_ume_001_step2
done

###############################################################################

echo
echo "10/14. Prep for step #3"

rm -rf $d0/step3

mkdir $d0/step3

cd $d0/step3

###############################################################################

echo
echo "11/14. Merge magn_orig (echo combined, coil uncombined)"

fslmerge -t magn_orig $d0/ch*/magn.nii*

###############################################################################

echo "12/Merge phase_orig (echo combined, coil uncombined)"

fslmerge -t phase_orig $d0/ch*/echofit_phase.nii*

###############################################################################

echo
echo 13/14. Create/import parameter files

echo 3.14159265359 > ptb_scaling_factor.txt

echo $deltaTE > ptb_TE.txt

cd $d0

if [ -f ch01/ptb_res_new.txt ]; then
	cp ch01/ptb_res_new.txt step3/ptb_res.txt
fi
cp -R dicom step3/
cp ptbs_ume_001_step3.m step3/
cp brainmask/roi.nii step3/
cp jmatcml step3/

###############################################################################

echo
echo "14/14. Coil combination, background removal and QSM inversion (step #3)"

cd $d0/step3

./jmatcml ptbs_ume_001_step3

###############################################################################

echo
echo End of script.
