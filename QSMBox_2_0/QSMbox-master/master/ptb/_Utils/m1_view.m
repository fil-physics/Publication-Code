function m1_view(froot,clipping_range,slice_range,panels)
% Display axial cuts
% 
% SINTAX
%  m1_view(fname)
%   Display NIFTI(.gz) image
% 
%  m1_view(fname,{cmin,cmax})
%   Clipping range options:  
%    ['gmin'] global min
%    ['gmax'] global max
%    ['gma2'] global max/2
%    ['mstd'] nonzero mean+-sigma (default)
% 
%  m1_view(fname,{cmin,cmax},slice_range)
%   Slice range:
%    [min,max], e.g. [0.2,0.75] (default) 
% 
%  m1_view(fname,{cmin,cmax},slice_range,panels)
%   Number of panels in montage:
%    [subplot_row,subplot_col], e.g. [3,6] (default)
%
% Created by Julio Acosta-Cabronero

if nargin < 2
    clipping_range = {'mstd','mstd'};
end

if nargin < 3
    slice_range = [.2,.75];
end

if nargin < 4
    panels = [3,6];
end

[mat] = m1_nifti_load(froot);

m1_fig(mat,clipping_range,slice_range,panels)
