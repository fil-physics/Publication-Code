function m1_fig(vol,clipping_range,slice_range,panels)
% DESCRIPTION
%  Display montage of 3D image matrix
%                   
% SINTAX
%  m1_fig(vol,{clipping_range},slice_range,panels)
%
% EXAMPLE
%  m1_fig(ima_mat,1.5,[.2,.75],[3,6])
%
% ARGUMENTS
%  vol                     Input 3D matrix
%  clipping_range          Signal intensity clipping range
%                            NOTE: cell input
%                            [1] {1x2 vector} -> {cmin, cmax}
%                            [2] {scalar}     -> <non0>(+-)scalar*sigma
%                               OTHER OPTIONS: 'gmin' (global min)
%                                              'gmax' (global max)
%                                              'gma2' (global max/2)
%                                   DEFAULT => 'mstd' (<non0>+-1.5*sigma)
%  slice_range             Slice range (percentage) to display
%                          [bottom_slice, top_slice]
%                           DEFAULT => [.2,.75] (20%-75% total # of
%                                                slices)
%                                      [0,1]    (display all)
%  panels                  Number of cuts
%                          [subplot_row, subplot_col]
%                           DEFAULT => [3,6] (3x6)
%
% Created by Julio Acosta-Cabronero

warning off

%% Preamble
% NaN to 0
vol(isnan(vol)==1)=0;
% Stats
gmin = min( vol(:) );
gmax = max( vol(:) );
thr = gmax/10;
gmean = mean( vol( vol>thr ));
gstd = std( vol( vol>thr ));
gmean(isnan(gmean)==1)=0;
gstd(isnan(gstd)==1)=0;

plane = 3;
matrix = size( vol );

%% Complete input arguments & defaults
% clipping range
if nargin < 2
    scalar = 1.5;
end

if exist( 'clipping_range' )
    if length( clipping_range ) == 1
        scalar = {clipping_range};
        clear clipping_range
    else
        scalar = 1.5;
    end
end

if ~exist( 'clipping_range' )
    clipping_range = { gmean-scalar*gstd, gmean+scalar*gstd };
end

cmin = cell2mat( clipping_range( 1 ));
cmax = cell2mat( clipping_range( 2 ));

if cmin == 'gmin'
    cmin = gmin;
end

if cmax == 'gmax'
    cmax = gmax;
end

if cmax == 'gma2'
    cmax = gmax/2;
end

if cmin == 'mstd'
    cmin = gmean-scalar*gstd;
end

if cmax == 'mstd'
    cmax = gmean+scalar*gstd;
end

if cmin == 0 && cmax == 0
    cmax = 1;
elseif cmin == 1 && cmax == 1
    cmin = 0;
end

% slice range
if nargin < 3
    slice_range = [ .2 .75 ];
end

if length( slice_range ) >= 2
    if max( slice_range(1:2) ) <= 1 && min( slice_range(1:2) ) >= 0 
        bottom_slice = slice_range( 1 );
        top_slice    = slice_range( 2 );
    else
        error('Incorrect slice range: input range [0 1]')
    end
else
    error('Incorrect argument: slice_range. See help.')
end

% panels
if nargin < 4
    panels = [ 3 6 ];
end

if length( panels ) >= 2
    subplot_row = panels( 1 );
    subplot_col = panels( 2 );    
else
    error('Incorrect argument: panels. See help.')
end
    
%% Crop volume
lb = round( bottom_slice*matrix( plane ));
if lb==0; lb = 1; end
ub = round( top_slice*matrix( plane ));
vol_crop = vol( :, :, lb:ub );

%% Select slices
matrix = size( vol_crop );
subplot_total = subplot_row*subplot_col;
slices = matrix( plane );
slices_to_display = round( slices * [ 1:subplot_total ]...
                                  * ( 1/( subplot_total+1 )));

%% Rotate volume
[ vol_disp ] = m1_rot_mat_fig( vol_crop );

%% Reshape volume
if plane == 3
    rsh_vol = reshape( vol_disp, [ matrix(2) matrix(1) 1 matrix(3)] );
else
    error( 'Unsupported display orientation' )
end
    
%% Display montage
montage( rsh_vol, 'Indices',        slices_to_display,...
                  'Size',           [subplot_row subplot_col],...
                  'DisplayRange',   [cmin cmax] )
impixelinfo
colorbar
drawnow

%% Print Stats
disp(' ')
disp([ 'IMAGE STATS' ])
disp([ ' global min   ' num2str(gmin) ])
disp([ ' global max   ' num2str(gmax) ])
disp([ ' nonzero mean ' num2str(gmean) ])
disp([ ' nonzero std  ' num2str(gstd) ])
disp([ ' ' ])
