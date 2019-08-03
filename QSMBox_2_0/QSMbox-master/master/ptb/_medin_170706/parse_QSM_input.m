function [lam RDF N_std Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir merit smv radius data_weighting gradient_weighting Debug_Mode] = parse_QSM_input(varargin)
% Internal function to parse input arguments (adapted from MEDI toolbox)
%   Created by Tian Liu on 2011.03.16
%   Modified by Tian Liu and Shuai Wang on 2011.03.15
%   Modified by Tian Liu on 2013.07.24
%   Last modified by Julio Acosta-Cabronero, 20 Jul 2017

merit               = 1;
smv                 = 0;
lam                 = 1000;
radius              = 5;
data_weighting      = 1;
gradient_weighting  = 1;
pad                 = 0;
matrix_size0        = 0;
Debug_Mode          = 'NoDebug';

filename = ['RDF.mat'];
if size(varargin,2)>0
    for k=1:size(varargin,2)
        if strcmpi(varargin{k},'filename')
            filename = varargin{k+1};
        end
        if strcmpi(varargin{k},'lambda')
            lam = varargin{k+1};
        end
        if strcmpi(varargin{k},'data_weighting')
            data_weighting = varargin{k+1};
        end
        if strcmpi(varargin{k},'gradient_weighting')
            gradient_weighting = varargin{k+1};
        end
        if strcmpi(varargin{k},'merit')
            merit = varargin{k+1};
        end
        if strcmpi(varargin{k},'smv')
            smv = varargin{k+1};
        end
        if strcmpi(varargin{k},'radius')
            radius = varargin{k+1};
        end
        if strcmpi(varargin{k},'pad')
            pad = varargin{k+1};
        end
        if strcmpi(varargin{k},'DEBUG')
            Debug_Mode = varargin{k+1};
        end
    end
end

load(filename,'RDF', 'N_std', 'Mask', 'matrix_size', 'voxel_size', 'delta_TE' ,'CF', 'B0_dir');

if exist('delta_TE','var')==0
    delta_TE = input('TE spacing = ');
    save(filename,'delta_TE','-append');
end

if exist('CF','var')==0
    CF = input('Central Frequency = ');
    save(filename,'CF','-append');
end

if exist('B0_dir','var')==0
    B0_dir = input('B0 direction = ');
    save(filename,'B0_dir','-append');
end

if sum(pad(:))
    matrix_size0 = matrix_size;
    matrix_size = matrix_size + pad;
    RDF = padarray(RDF, pad,'post');
    N_std = padarray(N_std, pad,'post');
    Mask = padarray(Mask, pad,'post');
end
end
