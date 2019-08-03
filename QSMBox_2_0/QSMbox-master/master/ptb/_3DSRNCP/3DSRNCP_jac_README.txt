  The 3DSRNCP unwrapper was originally written by the authors in C++.

  The code was modified by HF Sun to output 3D second differences.

  If you hit issues with the existing versions of 3DSRNCP, you can
  generate a binary executable for your own OS using Intel's C++ compiler.

  The C++ code can be found in the present directory.

  Install Intel C++ compiler from Parallel Studio XE:

   https://software.intel.com/en-us/qualify-for-free-software

  Compile, e.g.

   icpc -std=c++11 3DSRNCP_jac.cpp -o 3DSRNCP_newOS


  If you do this, you will also need to do the following:

  0. In the terminal window:

        cp 3DSRNCP_deb.m 3DSRNCP_newOS.m

  1. Add a new entry in 'm1_shf_3DSRNCP.m' ->[switch OS_type... end]

        switch OS_type
            case 2.1
                unwrapper = '3DSRNCP_deb';
            case 2.11
                unwrapper = '3DSRNCP_ubuntu-14.04';
            case 2.12
                unwrapper = '3DSRNCP_deb7_wheezy';
            case 2.2
                unwrapper = '3DSRNCP_rh';
            case 2.3
                unwrapper = '3DSRNCP_newOS';
        end

  2. Similarly, update 'QSMbox/master/ptb/_PhaseTools/m1_shf_mgre_comb.m'

        switch OS_type
            case 2.1
                unwrapper = '3DSRNCP_deb';
            case 2.11
                unwrapper = '3DSRNCP_ubuntu-14.04';
            case 2.12
                unwrapper = '3DSRNCP_deb7_wheezy';
            case 2.2
                unwrapper = '3DSRNCP_rh';
            case 2.3
                unwrapper = '3DSRNCP_newOS';
        end

  3. Update '~/.QSMbox/ptb_OS.m'

        function OS = ptb_OS
        OS = 2.30;

  4. Please share the binary(ies) and code via https://gitlab.com/acostaj/QSMbox
     
==============================================================================

Update from Renat Yakupov and Arturo Cardenas-Blanco, Magdeburg (27 Sep 2017):

libstdc++.so.6 and libc.so.6 too old

While running the QSMbox we ran into a problem of library incompatibilities with the following error messages:

libstdc++.so.6: version `GLIBCXX_3.4.21' not found
libc.so.6: version `GLIBC_2.14' not found

A solution suggested in QSMbox/master/ptb/_3DSRNCP/3DSRNCP_jac_README.txt didnt work as we dont have the rights to implement it. Instead, Arturo found the following workaround:

I. To get the latest version of libc.so.6, perform the following steps:

1. mkdir ~/glibc_install; cd ~/glibc_install
2. wget http://ftp.gnu.org/gnu/glibc/glibc-2.14.tar.gz
3. tar zxvf glibc-2.14.tar.gz
4. cd glibc-2.14
5. mkdir build
6. cd build
6. ../configure --prefix=/home/<USERNAME>/glibc-2.14
8. make -j 4
9. make install

II. To get the latest version of libstdc++.so.6, install the latest version of Miniconda from https://conda.io/miniconda.html.

III. Add the following line to Matlab startup file:

setenv('LD_LIBRARY_PATH', '/home/<USERNAME>/miniconda3/lib');

IV. Create a link to libc.so.6 in miniconda3/lib directory:

ln -s /home/<USERNAME>/glibc-2.14/lib/libc.so.6 /home/<USERNAME>/miniconda3/lib/

