#!/bin/bash
#++++++++++++++++++++++++++ Script begins here ++++++++++++++++++++++++++
if [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ]; then
        ########## HELP begins here ##########
        echo
        echo "Run MATLAB program from shell"
        echo " \$ ./`basename $0` <script name>"
	echo
	echo "Copy it to e.g. /usr/local/bin and run it as:"
	echo " \$ `basename $0` <script name>"
	echo
	echo "SCRIPT NAME"
	echo " fileroot, no extension"
	echo " script must be in pwd"
        echo
	echo Created by Julio Acosta-Cabronero
        ########### HELP ends here ###########
else
        if [ -z "$1" ]; then
                $0 --help
        else
                ########## BODY begins here ##########
                if [ -n "$1" ]; then
			# Specify MATLAB version (edit below if required):
                          MATLABCALL=matlab

                        # Record date ID
			  TS=`date +"%Y%m%d_%H%M%S"`

			# Stamp date ID on log’s filename and print full path
                          echo "To monitor progress, open new terminal window (Ctrl+Shift+t) and paste:"
                          echo "tail -f $PWD/matop_"$TS"_"$1".log"

			# Run MATLAB silently 
                          time nohup ${MATLABCALL} -nodesktop -nosplash -r "${1}; quit" < /dev/null > matop_"$TS"_"$1".log 2>&1
                fi
                ########### BODY ends here ###########
        fi
fi
#+++++++++++++++++++++++++++ Script ends here +++++++++++++++++++++++++++
