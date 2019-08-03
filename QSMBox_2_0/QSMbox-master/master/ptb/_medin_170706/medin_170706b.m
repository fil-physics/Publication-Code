function [x,cost_reg,cost_fid,N_CG_iter,CG_residual] = medin_170706b(varargin)
% Non-linear morphology enabled dipole inversion (nMEDI)
%   [x,cost_reg,cost_fid] = medin_170706b(varargin)
%
%   OUTPUT
%   x               susceptibility distribution 
%   cost_reg        cost of the regularization term
%   cost_fid        cost of the data fidelity term
%   N_CG_iter       # of CG iterations per qN iteration
%   CG_residual     CG residuals across qN iters
%   
%   INPUT
%   RDF.mat in pwd  
%   medin_170706b('lambda',lam,...)      lam specifies the regularization parameter
%                                         lam weights the data fidelity term
%
%   ----optional----   
%   medin_170706b('smv',[0/1],...)       SMV deconvolution (disable/enable)
%   medin_170706b('radius',radius,...)   specify the radius for the SMV
%                                         operator using differential form
%   medin_170706b('merit',[0/1],...)     model error reduction through
%                                         iterative tuning (disable/enable)
%   medin_170706b('pad',padsize,...)     zeropad
%
%   REFERENCES 
%   T. Liu et al. MRM 2013;69(2):467-76
%   J. Liu et al. Neuroimage 2012;59(3):2560-8.
%   T. Liu et al. MRM 2011;66(3):777-83
%   de Rochefort et al. MRM 2010;63(1):194-206
%
%   Adapted from Cornell's MEDI toolbox
%   Last modified by Julio Acosta-Cabronero

% RDF.mat: RDF N_std Mask matrix_size voxel_size delta_TE CF B0_dir 
[lambda,RDF,N_std,Mask,matrix_size,matrix_size0,voxel_size,...
 delta_TE,CF,B0_dir,merit,smv] = parse_QSM_input(varargin{:});

% genDI.mat: stp ptb e cg_max_iter cg_tol max_iter tol_norm_ratio
%            data_weighting_mode gradient_weighting_mode wC SMV_radius
%            merit_f merit_p MSDI_flag downsampl_factor gadmpm_flag
load genDI

disp(['Lambda (consistency multiplier) = ' num2str(round(lambda))])

% Operator definitions
grad = @m1_tl_medin_fgrad;
div  = @m1_tl_medin_bdiv;

% Dipole kernel
disp(['Calculate dipole kernel for B0_dir = [' num2str(B0_dir) ']'])
D = m1_tl_medin_dipole_kernel(matrix_size,voxel_size,B0_dir);

if smv
    disp('Enable SMV deconvolution')
    S = m1_tl_medin_SMV_kernel(matrix_size,voxel_size,SMV_radius);
    D = S.*D;
end

% Initialise tempn as phase sigma
if data_weighting_mode==0; N_std = 1; end
tempn = double(N_std).*Mask;

% Consistency weighting
m = m1_tl_medin_dataterm_mask(data_weighting_mode,tempn,Mask); % transparent if data_weighting_mode=0;
% if data_weighting_mode==1
disp('Compute consistency weighting')
m = m.*wC;
if merit; data_weighting_mode = 1; end 
N_std = Mask./m; N_std(isnan(N_std))=0; N_std(isinf(N_std))=0; 
% end
tempn = double(N_std).*Mask;
% ttt = m; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
%  m1_nifti_save(ttt,stp,'W_fid'); clear ttt
% end

% Regulariser spatial constrain
switch gradient_weighting_mode
case 0
    disp('Skip grad mask')
    wG = 1;
     
case 1
     [ima] = m1_nifti_load('edgeima');
    if gadmpm_flag; ttt = ima; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); ima = ttt; clear ttt; end
    if downsampl_factor>1
        ima = m1_downsampl(ima,downsampl_factor);
    end
    if matrix_size0
        ima = padarray(ima,matrix_size-matrix_size0,'post');
    end
    disp('Compute gradient mask')
    wG = m1_tl_medin_gradient_mask(gradient_weighting_mode,ima,Mask,grad,voxel_size,0.9); clear ima
    wGx=wG(:,:,:,1);               wGy=wG(:,:,:,2);              wGz=wG(:,:,:,3);
    SavewGxyz = false;
    if SavewGxyz
        ttt = wGx; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'wGx');
        ttt = wGy; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGy'); 
        ttt = wGz; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGz'); 
        ttt = wGx+wGy+wGz; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGsum')
    end
     clear wGx wGy wGz ttt

case 2
     [ima] = m1_nifti_load('edgeima');
    if gadmpm_flag; ttt = ima; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); ima = ttt; clear ttt; end
    if downsampl_factor>1
        ima = m1_downsampl(ima,downsampl_factor);
    end
    if matrix_size0
        ima = padarray(ima,matrix_size-matrix_size0,'post'); 
    end
    disp('Compute gradient mask')
    wG = m1_tl_medin_gradient_mask(gradient_weighting_mode,ima,Mask,grad,voxel_size,0.9); clear ima
    wGx=wG(:,:,:,1);               wGy=wG(:,:,:,2);              wGz=wG(:,:,:,3);
    if data_weighting_mode==1 && ( length(wC(:))>1 && min(wC(Mask>0))<1 )
        disp('Force regularisation of unreliable phases')
        wGx = wGx+(1-wC); wGx(wGx>0)=1;
        wGy = wGy+(1-wC); wGy(wGy>0)=1;
        wGz = wGz+(1-wC); wGz(wGz>0)=1;
    end
    SavewGxyz = false;
    if SavewGxyz
        ttt = wGx; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'wGx');
        ttt = wGy; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGy'); 
        ttt = wGz; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGz'); 
        ttt = wGx+wGy+wGz; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end 
         m1_nifti_save(ttt,stp,'wGsum'); clear ttt
    end
    wG(:,:,:,1)=wGx;               wG(:,:,:,2)=wGy;              wG(:,:,:,3)=wGz;
    clear wGx wGy wGz
end

clear wC

% MERIT on?
if merit
    disp('MERIT enabled')
end

% Complex data (measured)
b0 = m.*exp(1i*RDF);

% Initialise iterative loop
iter                    = 0;
x                       = zeros(matrix_size);
res_norm_ratio          = Inf;
res_norm_ratio_history  = zeros(1,max_iter);
cost_fid_history        = zeros(1,max_iter);
cost_reg_history        = zeros(1,max_iter);
cost_reg                = 0;
cost_fid                = 0;

while res_norm_ratio>tol_norm_ratio && iter<max_iter
        
    iter=iter+1;
    
    % conjugate gradient computation
    Ln = 1;
    if Ln==1
        Vr      = 1./sqrt(abs(wG.*grad(real(x),voxel_size)).^2+e);
    elseif Ln==2
        Vr      = 2;
    end
    w           = m.*exp(1i*ifftn(D.*fftn(x)));
    reg         = @(dx) div(wG.*(Vr.*(wG.*grad(real(dx),voxel_size))),voxel_size);
    fidelity    = @(dx) 2*lambda*real(ifftn(D.*fftn(conj(w).*w.*real(ifftn(D.*fftn(dx))))));
    
    A           = @(dx) reg(dx)+fidelity(dx);
    b           = reg(x) + 2*lambda*real(ifftn(D.*fftn(conj(w).*conj(1i).*(w-b0))));
    
    [dx,cg_rsd,N_iter]   = m1_tl_medin_cgsolve(A,-b,cg_tol,cg_max_iter,0);
    
    dx                  = real(dx);
    N_CG_iter(iter)     = N_iter;
    CG_residual(iter)   = cg_rsd;
    
    % relative update size 
    res_norm_ratio = norm(dx(:))/norm(x(:));
        
    % update solution
    x           = x+dx;
    
    % save intermediate solution
    if iter==1
        SaveIter    = true;
        default_tnr = 0.1;
    end
    if SaveIter && res_norm_ratio<default_tnr && tol_norm_ratio<default_tnr
        SaveIter = false;
        x0 = x/(2*pi*delta_TE*CF)*1e6.*Mask;
        if matrix_size0
            x0 = x0(1:matrix_size0(1),1:matrix_size0(2),1:matrix_size0(3));
        end  
        SaveHistory = false;
        if SaveHistory
            ttt = x0; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,['x_iter' num2str(iter)]);
        else
            ttt = x0; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'x');
        end
        clear x0 ttt
    end
        
    % residual matrix
    wres        = m.*exp(1i*(real(ifftn(D.*fftn(x))))) - b0;

    % calculate costs & archive history
    cost_fid                        = norm(wres(:),2);
    cost_fid_history(iter)          = cost_fid;
    cost_reg                        = abs(wG.*grad(x));
    cost_reg                        = sum(cost_reg(:));
    cost_reg_history(iter)          = cost_reg;
    res_norm_ratio_history(iter)    = res_norm_ratio;

    % MERIT
    if merit
        wres            = wres-mean(wres(Mask(:)==1));
        a               = wres(Mask(:)==1);
        factor          = std(abs(a))*merit_f;
        wres            = abs(wres)/factor;
        wres(wres<1)    = 1;
        N_std(Mask==1)  = N_std(Mask==1).*wres(Mask==1).^merit_p;
        tempn           = N_std;
        m               = m1_tl_medin_dataterm_mask(data_weighting_mode,tempn,Mask);
        b0              = m.*exp(1i*RDF);
        SaveMERITtemp = false;
        if SaveMERITtemp
            ttt = m; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'W_merit'); clear ttt
        end
    end
    
    disp(' ')
    verbose = 1;
    if smv && verbose
        disp(['++Differential SMV (kernel radius=' num2str(SMV_radius) ' mm)'])
    end
    fprintf('++%d.  Update = %.3f (tol=%.0e).  Consistency = %.2e.  Regulariser = %.2e\n',...
             iter,res_norm_ratio,tol_norm_ratio,cost_fid,cost_reg);
    disp(' ')
    toc %
    
%     if MSDI_flag==-1
%         if N_iter==cg_max_iter
%             break
%         end
%     end    
end

% save MERIT adjustment to consistency weighting
SaveMERIT = false;
if SaveMERIT
    ttt = m; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'W_merit'); clear ttt
end

% save history to text
SaveHistory = false;
if SaveHistory
    save('update_history.txt',  'res_norm_ratio_history','-ascii')
    save('cost_fid_history.txt','cost_fid_history',      '-ascii')
    save('cost_reg_history.txt','cost_reg_history',      '-ascii')
end

% convert x to ppm
x = x/(2*pi*delta_TE*CF)*1e6;

% crop
if matrix_size0
    x = x(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
end

end
