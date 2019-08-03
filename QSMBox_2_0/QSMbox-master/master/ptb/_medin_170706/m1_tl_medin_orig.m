function [ x_froot, cost_fid, cost_reg, N_brain ] =...
            m1_tl_medin( B_local_froot, W_froot, lambda,...
                       data_weighting, gradient_weighting,...
                       merit, smv, radius, pad,...
                       cg_max_iter, cg_tol, max_iter, tol_norm_ratio, ...
                       res, B0, TE )
% DESCRIPTION
%  Non-linear, L1-norm regularised, 
%   morphology-enabled magnetic dipole inversion (MEDIN)
% 
% SINTAX
%  [ x_froot, cost_fidelity, cost_regularisation, brainmask_size ] =...
%             m1_tl_medin( B_local_froot, W_froot, lambda,...
%                        data_weighting, gradient_weighting,...
%                        MERIT, SMV, radius, padsize,...
%                        cg_max_iter, cg_tol, max_iter, tol_norm_ratio, ...
%                        res, B0, TE )
%
% OUTPUT
%  x - QSM as NIFTI(.gz)
%  x_froot - QSM fileroot
%  cost_fid - the cost of the regularization term
%  cost_reg - the cost of the data fidelity term
%  N_brain - number of nonzero voxels in weighting matrix
%
% When using the code, please cite:
%  T. Liu et al. MRM 2013;69(2):467-76
%  J. Liu et al. Neuroimage 2012;59(3):2560-8.
%  T. Liu et al. MRM 2011;66(3):777-83
%  de Rochefort et al. MRM 2010;63(1):194-206
%
% > Adapted from Ildar Khalidov
% >> Modified by Tian Liu on 2011.02.01
% >>> Modified by Tian Liu and Shuai Wang on 2011.03.15
% >>>> Modified by Tian Liu and Shuai Wang on 2011.03.28
% >>>>> Modified by Tian Liu on 2013.07.24
% >>>>>> Last modified by Julio Acosta-Cabronero on 2014.01.20

%% Complete input
c = 1;
if nargin < c
    [ fname, B_local_path ] = uigetfile(...
                                { '*.nii.gz'; '*.nii' },...
                                'Select local induction (B_local) image' );
    [ B_local_froot ] = m1_fname_split( fname, 'gz' );
end

if exist( 'B_local_path' )
    cd( B_local_path )
end

c = c+1;
if nargin < c
    [ fname ] = uigetfile(...
                                   { '*.nii.gz'; '*.nii' },...
                                   'Select noise weighting image' );
    [ W_froot ] = m1_fname_split( fname, 'gz' );
end

c = c+1;
if nargin < c,
    lambda = 1000;
end

c = c+1;
if nargin < c
    data_weighting = 1;
end

c = c+1;
if nargin < c
    gradient_weighting = 1;
end

c = c+1;
if nargin < c
    merit = 1;
end

c = c+1;
if nargin < c
    smv = 1;
end

c = c+1;
if nargin < c
    radius = 4;
end

c = c+1;
if nargin < c
    pad = 0; %[100 100 50];
end

c = c+1;
if nargin < c
    cg_max_iter = 1000;
end

c = c+1;
if nargin < c
    cg_tol = 0.01;
end

c = c+1;
if nargin < c
    max_iter = 50;
end

c = c+1;
if nargin < c
    tol_norm_ratio = 0.1;
end

c = c+1;
if nargin < c
    if exist( 'dicom' )
        [ res, B0, TE ] = m1_dicom_extract_res_B0_TE( 'dicom' );
    end
    
    if exist( 'res.mat' )
        ttt = load('res.mat'); res = ttt.res; clear ttt
    end
    if exist( 'B0.mat' ); load B0; end
    if exist( 'TE.mat' ); load TE; end
    
    if ~exist( 'res' ); if exist( 'res.txt' )
        res = double(load('res.txt'));
    else
        res = input( '> Voxel resolution (vector in mm) => e.g. [1 1 2]: ' );
    end; end
    if ~exist( 'B0' ); if exist( 'B0.txt' )
        B0 = double(load('B0.txt'));
    else
        B0 = input( '> Field strength (in T): ' );
    end; end
    if ~exist( 'TE' ); if exist( 'TE.txt' )
        TE = double(load('TE.txt'));
    else
        TE = input( '> Echo time (in s): ' );
    end; end
    if ~exist( 'B0_dir' ); if exist( 'B0_dir.txt' )
        B0_dir = double(load('B0_dir.txt'));
    else
        B0_dir = input( '> B0 direction (Z cosines) => e.g. [0 0 1] for true axial: ' );
    end; end
end

viewdisp=0;
saveall=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('> Non-linear MEDI')
Debug_Mode = 'NoDebug';

%% Compute constants
gyro = 2*pi*42.576;
ppm2rad = double(B0*gyro*TE(1));
rad2ppm = double(1/ppm2rad);

save('ppm2rad.txt','ppm2rad','-ascii')
save('rad2ppm.txt','rad2ppm','-ascii')

res
B0
TE
ppm2rad
rad2ppm

%% Load data
[ RDF, stp ] = m1_nifti_load( B_local_froot );
[ iMag ] = m1_nifti_load( W_froot );

%% Prep data
RDF = ppm2rad*RDF;
mask = iMag; mask(iMag>0) = 1;
RDF = RDF.*mask;
N_brain = numel( mask(iMag>0) );
N_std = 1./iMag;
N_std(isnan(N_std)) = 0;
N_std(isinf(N_std)) = 0;
N_std = N_std .* mask;
matrix_size = size( RDF );

%% Zeropad
if sum(pad(:))
    RDF = padarray(RDF, pad);
    N_std = padarray(N_std, pad);
    iMag = padarray(iMag, pad);
    mask = padarray(mask, pad);
    matrix_size = size( RDF );
end

%% Weights definition 
data_weighting_mode = data_weighting;
gradient_weighting_mode = gradient_weighting;
%grad = @m1_tl_medin_cgrad;
%div = @m1_tl_medin_cdiv;
grad = @m1_tl_medin_fgrad;
div = @m1_tl_medin_bdiv;

%% Summary dir
DIR = [ 'MEDIN_progress/lambda' num2str(lambda) '_'...
                       'tnr' num2str(tol_norm_ratio) '_'...
                       'cg_iter' num2str(cg_max_iter) '_'...
                       'cg_tol' num2str(cg_tol) '_'...
                       'MERIT' num2str(merit) '_'...
                       'SMV' num2str(smv) '_'...
                       'radius' num2str(radius) '_'...
                       'W' num2str(data_weighting_mode) '_'...
                       'Mg' num2str(gradient_weighting_mode) '_'...
                       B_local_froot ];

mkdir( DIR )
disp([ '> Make directory: ' DIR ])

%% Dipole kernel
disp('> Compute dipole kernel')
D = m1_tl_medin_dipole_kernel( matrix_size, res, B0_dir );

tempn = double(N_std);

%% SMV
if (smv)
    disp('> SMV')
    S = m1_tl_medin_SMV_kernel( matrix_size, res, radius );
    mask = m1_tl_medin_SMV( mask, matrix_size, res, radius )>0.999;
    D = S.*D;
    RDF = RDF - m1_tl_medin_SMV( RDF, matrix_size, res, radius );
    RDF = RDF .* mask;
    tempn = sqrt(...
          m1_tl_medin_SMV( tempn.^2, matrix_size, res, radius )+tempn.^2 );

    % display
    if viewdisp==1
    figure(918)
    m1_figure_display_mat(RDF,{-0.03*ppm2rad,0.03*ppm2rad})
    title('SMV-filtered RDF')    
    saveas(918, 'fig918_SMV_RDF', 'png')
    end
    
    % save data
    m1_nifti_save(rad2ppm*RDF,stp,['SMV' num2str(radius) '_' B_local_froot]);
end

%% Noise weighting mask
disp('> Dataterm mask')
m = m1_tl_medin_dataterm_mask( data_weighting_mode, tempn, mask );
b0 = m.*exp(1i*RDF);

    % display
    if viewdisp==1
    figure(919)
    m1_figure_display_mat(m)
    title('Dataterm mask (init)')    
    saveas(919, 'fig919_dataterm_init', 'png')
    end

%     figure(920)
%     m1_figure_display_mat(b0)
%     title('b0 (init)')    
%     saveas(920, 'fig920_b0_init', 'png')
   
    if saveall==1
    % save data
    m1_nifti_save(m,stp,'dataterm_init');
    m1_nifti_save(angle(b0),stp,'b0_init');
    end  

%% Gradient mask
disp('> Gradient mask')
Mg = m1_tl_medin_gradient_mask( gradient_weighting_mode,...
                                iMag, mask, grad, res);

    % display
    if viewdisp==1
    figure(921)
    m1_figure_display_mat(Mg,{0,1})
    title('Gradient mask')    
    saveas(921, 'fig920_gradient_mask', 'png')
    end

%     % save data
     %m1_nifti_save(squeeze(Mg(:,:,:,1)),stp,'Mgx');
     %m1_nifti_save(squeeze(Mg(:,:,:,2)),stp,'Mgy');
     %m1_nifti_save(squeeze(Mg(:,:,:,3)),stp,'Mgz');

%% Pre-loop
close all

x = zeros(matrix_size); %real(ifftn(conj(D).*fftn((abs(m).^2).*RDF)));
if (~isempty(findstr(upper(Debug_Mode),'SAVEITER')))
    m1_tl_medin_store_CG_results((rad2ppm*x).*mask);
end

if merit
    disp('> MERIT on')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CG LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

rnr = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);
e=eps;
badpoint = zeros(matrix_size);
iter=0;
% m1_sfigure(920)   
disp('> Initialise loop...')
while (rnr>tol_norm_ratio)&&(iter<max_iter)
    iter=iter+1;
    
    % conjugate gradient computation
    Vr = 1./sqrt(abs(Mg.*grad(real(x),res)).^2+e);
    w = m.*exp(1i*ifftn(D.*fftn(x)));
    
    reg = @(dx) div(Mg.*(Vr.*(Mg.*grad(real(dx),res))),res);
    fidelity = @(dx) 2*lambda*real(ifftn(D.*fftn(conj(w).*...
                                            w.*real(ifftn(D.*fftn(dx))))));

    A = @(dx) reg(dx)+fidelity(dx);       
    b = reg(x) + 2*lambda*real(ifftn(D.*fftn( conj(w).*conj(1i).*(w-b0))));

    dx = real( m1_tl_medin_cgsolve( A, -b, cg_tol, cg_max_iter, 0 ));
    
    % compute residual ratio
    rnr = norm(dx(:))/norm(x(:));
    res_norm_ratio(iter) = rnr;
    
    % update susceptibility distribution
    x = x + dx;
    
    % display progress
    if viewdisp==1
    m1_sfigure( 950 )
    m1_figure_display_mat( x, {-0.08*ppm2rad, 0.08*ppm2rad} )
    title([ 'Susceptibility Regularisation Progress (lambda=' num2str(lambda) ') - Iteration #' num2str(iter) ' (residual=' num2str(rnr) ')' ])
    eval([ 'saveas(950, ''fig950_QSM_' num2str(iter+10) ''', ''png'')' ])
    %saveas( 950, 'fig950_QSM_progress', 'png' )
    end

    %if saveall==1
    % save progress
    x_prog = (rad2ppm*x) .* mask;
    cd( DIR )
    %m1_nifti_save( x_prog, stp, 'QSM_progress' );
    m1_nifti_save( x_prog, stp, [ 'QSM_progress_' num2str(iter+10) ]);
    cd ../..
    %end

    % compute residual matrix
    wres=m.*exp(1i*(real(ifftn(D.*fftn(x))))) - b0;

    cost_data_history(iter) = norm( wres(:), 2 );
    cost = abs( Mg.*grad(x) );
    cost_reg_history(iter) = sum( cost(:) );
    
    % MERIT
    if merit
        wres = wres - mean(wres(mask(:)==1));
        a = wres(mask(:)==1);
        factor = std(abs(a))*6;
        wres = abs(wres)/factor;
        wres(wres<1) = 1;
        badpoint(wres>1)=1;
        N_std(mask==1) = N_std(mask==1).*wres(mask==1).^2;
        tempn = N_std;
        if (smv)
         tempn=sqrt(...
           m1_tl_medin_SMV(tempn.^2, matrix_size, res, radius) + tempn.^2);
        end
        m = m1_tl_medin_dataterm_mask(data_weighting_mode, tempn, mask);
        b0 = m.*exp(1i*RDF);
    end
    
    % print progress
    fprintf('=> iter: %d; res_norm_ratio: %8.4f; cost_L2: %8.4f; cost_reg: %8.4f.\n',iter,rnr,cost_data_history(iter),cost_reg_history(iter));
toc
end
disp('...CG loop terminated.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display & save CG loop final updates
% m
if viewdisp==1
figure(919)
m1_figure_display_mat(m)
title('Dataterm mask')    
saveas(919, 'fig919_dataterm_final', 'png')
end

if saveall==1
m1_nifti_save(m,stp,'dataterm_final');
% b0
% figure(920)
% m1_figure_display_mat(b0)
% title('b0 (final)')    
% saveas(920, 'fig920_b0_final', 'png')
m1_nifti_save(b0,stp,'b0_final');
end

%% Summary
cost_fid = cost_data_history(iter);
cost_reg = cost_reg_history(iter);

%% mask x & convert to ppm
x = (rad2ppm*x) .* mask;

%% Display montage
if viewdisp==1
ttt = '> Final QSM';
disp(ttt)
m1_sfigure( 999 )
m1_figure_display_mat( x, {-0.08,0.08} )
title(ttt)
saveas(999, 'fig999_QSM', 'png')
end

%% Crop (if padding enabled)
if sum(pad(:))
  x = x(pad(1)+1:end-pad(1),pad(2)+1:end-pad(2),pad(3)+1:end-pad(3));
  iMag = iMag(pad(1)+1:end-pad(1),pad(2)+1:end-pad(2),pad(3)+1:end-pad(3));
  RDF = RDF(pad(1)+1:end-pad(1),pad(2)+1:end-pad(2),pad(3)+1:end-pad(3));
  mask = mask(pad(1)+1:end-pad(1),pad(2)+1:end-pad(2),pad(3)+1:end-pad(3));
end

%% Save data
%m1_tl_medin_store_QSM_results(x, iMag, RDF, mask,...
m1_tl_medin_store_QSM_results('Norm', 'L1','Method','MEDIN','Lambda',lambda,...
                  'SMV',smv,'Radius',radius,'IRLS',merit,...
                  'res',res,'matrix_size',matrix_size,...
                  'Data_weighting_mode',data_weighting_mode,...
                  'Gradient_weighting_mode',gradient_weighting_mode,...  
                  'L1_tol_ratio',tol_norm_ratio, 'Niter',iter,...
                  'CG_tol',cg_tol,'CG_max_iter',cg_max_iter,...
                  'B0_dir', B0_dir);

lambda=double(lambda);
save([DIR '/res_norm_ratio.mat'], 'res_norm_ratio' )
save([DIR '/cost_data_history.mat'], 'cost_data_history' )
save([DIR '/cost_reg_history.mat'], 'cost_reg_history' )
save([DIR '/lambda.txt'], 'lambda', '-ascii' )

% save NIFTi(.gz)
x_froot = [ 'QSM_MEDIN_lambda' num2str(lambda) '_'...
                       'tnr' num2str(tol_norm_ratio) '_'...
                       'cg_iter' num2str(cg_max_iter) '_'...
                       'cg_tol' num2str(cg_tol) '_'...
                       'MERIT' num2str(merit) '_'...
                       'SMV' num2str(smv) '_'...
                       'radius' num2str(radius) '_'...
                       'W' num2str(data_weighting_mode) '_'...
                       'Mg' num2str(gradient_weighting_mode) '_'...
                       B_local_froot ];
                                                   
m1_nifti_save(x,stp,x_froot);

%pause(10)

end
