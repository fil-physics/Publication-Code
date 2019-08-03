function [ res, B0, TE ] = m1_dicom_extract_res_B0_TE( DicomFolder )
% DESCRIPTION
%  DICOM field extractor
% 
% SINTAX
%  [ res, B0, TE ] = m1_dicom_extract_res_B0_TE( DicomFolder )
%
% INPUT
%  DicomFolder: It should only contain DICOM files
%               Present outputs can be generated from a single DICOM
%
% OUTPUT
%  res: Voxel resolution in mm
%  B0:  Field strength in T
%  TE:  Echo time in s
%
% > created by Shuai Wang (19/7/2013)
% >> adapted by Tian Liu (Read_Siemens_DICOM in MEDI toolbox)
% >>> last modified by Julio Acosta-Cabronero (14/10/2016)

if ~isdir( DicomFolder )
	disp(['ERROR: DICOM dir, ' DicomFolder ', does not exist'])
	res=false; B0=false; TE=false;
	return
end

    filelist = dir( DicomFolder );
    x = 1;
    while x <= length( filelist )
        if filelist(x).isdir == 1
            filelist = filelist([ 1:x-1 x+1:end ]);   % skip folders
        else
            x = x+1;
        end
    end

    ttt = [ DicomFolder '/' filelist(1).name ];
    info = dicominfo( ttt );
    res(1,1) = single( info.PixelSpacing(1) );
    res(1,2) = single( info.PixelSpacing(2) );
    res(1,3) = single( info.SliceThickness );
    CF = info.ImagingFrequency*1e6;
    NumEcho = info.EchoNumber;
    TE = single( zeros([ NumEcho 1 ]));
    if TE( info.EchoNumber )==0
        TE( info.EchoNumber ) = info.EchoTime*1e-3;
    end
    TE=TE(TE>0);
    TE=TE(1);
    B0 = CF/(42.57747892*1e6);
   
    save( 'dicom_info.mat', 'info' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matrix_size(1) = single(info.Width);
    matrix_size(2) = single(info.Height);    
    minSlice = 1e10;
    maxSlice = -1e10;
        if info.SliceLocation<minSlice
            minSlice = info.SliceLocation;
            minLoc = info.ImagePositionPatient;
        end
        if info.SliceLocation>maxSlice
            maxSlice = info.SliceLocation;
            maxLoc = info.ImagePositionPatient;
        end
    matrix_size(3) = round(norm(maxLoc - minLoc)/res(3)) + 1;    
    Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*res(3))];
    if matrix_size(3) ~= 1
        B0_dir = Affine3D\[0 0 1]';
    else
        if Affine2D\[0 0 1]'==[0;0]
            B0_dir = [0;0;1];
        else
            a=Affine2D\[0 0 1]';
            z=sqrt(1-a(1)^2-a(2)^2);
            B0_dir = [a;z];
        end
    end
    if B0_dir~=[0;0;1]
        disp('> Warning: Third dimension not perpendicular to B0!')
        !touch WARNING__SLICES_NOT_STRAIGHT_AXIAL
        B0_dir;
    else
        B0_dir;
    end

    B0_dir=double(abs(B0_dir)');
    if ~exist('ptb_B0_dir.txt')
        save('ptb_B0_dir.txt','B0_dir','-ascii')
    end
    CF=double(CF);
%     save('ptb_CF.txt','CF','-ascii')    
    B0=double(B0);
%     save('ptb_B0.txt','B0','-ascii')
    TE=double(TE);
%     save('ptb_TE.txt','TE','-ascii')
    res=double(res);
%     save('ptb_res.txt','res','-ascii')

end
