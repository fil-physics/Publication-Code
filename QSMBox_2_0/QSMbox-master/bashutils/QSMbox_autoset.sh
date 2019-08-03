#!/bin/bash
BASHPROF=${HOME}/.bash_profile

echo "Autoset QSMbox environment on login/terminal init ($BASHPROF)"
if [ -f ${BASHPROF} ]; then
    FLAG1=`cat ${BASHPROF} | grep mpb_profile.sh`
else
    FLAG1=
fi
if [ -z "$FLAG1" ]; then
        echo Update ${BASHPROF}
	echo " " >> ${BASHPROF}
        echo "# QSMbox" >> ${BASHPROF}
        echo "source ${PWD}/QSMbox_profile.sh" >> ${BASHPROF}
	source ${BASHPROF}
else
        echo QSMbox entry exists - do nothing
fi
