function [froot,ext] = m1_fname_split(fname,gz_flag)
% SINTAX
%  [fileroot,extension] = m1_fname_split(filename,'gz')
%
% OPTIONAL ARGUMENT
%  'gz': ignore '.gz' extension
%
% NOTES
%  This program searches for the last '.' in the filename
%  If several dots are present they will be included in the output fileroot
%  '.gz', however, can be ignored as a file extension using the argin 'gz' 
%
% Created by Julio Acosta-Cabronero

idx     = find(fname=='.');
froot   = fname(1:idx(end)-1);
ext     = fname(idx(end)+1:end );

if nargin==2
    if strcmpi(gz_flag,'gz') == 1
        if strcmpi(ext,'gz') == 1
            [froot,ext] = m1_fname_split(froot);
        end
    else
        error('Second argument must be ''gz'' ')
        return
    end    
end
