function m1_gzip(file_extension)
% SINTAX
%  m1_gzip(file_extension)
%
% Created by Julio Acosta-Cabronero

if nargin < 1
    file_extension = 'nii';
end

eval(['flag=dir(''*.' file_extension ''');']);

if ~isempty(flag)
    disp('Gzip')
    eval(['gzip(''*.' file_extension ''')']);
    eval(['delete(''*' file_extension ''')']);
else
    disp('Did nothing. Exit.')
end
