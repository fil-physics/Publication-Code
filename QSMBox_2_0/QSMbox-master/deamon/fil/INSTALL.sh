#!/bin/bash

MATSTART=/root/matlab/startup.m

d0=$(dirname `cd \`dirname "$0"\`; pwd`/`basename "$0"`)

cd $d0/../../..

d1=$PWD

echo
echo Edit $MATSTART
echo >> $MATSTART
echo "%% QSMbox" >> $MATSTART
echo addpath $d1/QSMbox >> $MATSTART
echo addpath $d1/QSMbox/master >> $MATSTART

echo
echo Update QSMbox/bashutils/QSMbox_local_settings.sh
echo Set variables: OS, SPMPATH
cp -v $d0/QSMbox_local_settings.sh $d1/QSMbox/bashutils/QSMbox_local_settings.sh

echo
echo Execute ./QSMbox/bashutils/QSMbox_local_settings.sh
cd $d1/QSMbox/bashutils
chmod +x QSMbox_local_settings.sh
./QSMbox_local_settings.sh

echo
echo Clone bashmritools
cd $d1
git clone https://gitlab.com/acostaj/bashmritools.git

echo 
echo Install pigz
apt-get install pigz

echo
echo Update bashmritools/src/gsettings
cp -v $d0/bashmritools_gsettings $d1/bashmritools/src/gsettings

echo
echo Execute ./bashmritools/INSTALL.sh
cd $d1/bashmritools
chmod +x INSTALL.sh
./INSTALL.sh

echo
echo Finished.
