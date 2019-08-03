function [index_opt,lambda_opt,Kappa] = m1_LcurveOpt(consistency, regularization, tauList, line_color)
% Written by: Job Bouwman,
% (modified code of Berkin Bilgic)
% contact: jgbouwman@hotmail.com
% date: 21-10-2014
% 
% cubic spline differentiation to find Kappa (largest curvature) 
% Selects the optimal parameter value for regularization form an L-curve.
% This should lead to an optimal balance in the reconstruction between
% data consistency and smoothness of solution.
% 
% See also: Hansen et al, SIAM J. Sci. Comput., 14(6), 1487–1503.
% "The Use of the L-Curve in the Regularization of Discrete Ill-Posed 
% Problems"
% 
% Last modified by Julio Acosta-Cabronero

% show the L-curve?   
DISPLAY_L_CURVE = true;  % line_color  = 'b';

% Detecting the valid nodes:
validIndices = and((consistency>0),(regularization>0));

% plotting (valid nodes in blue, invalid in red):
if DISPLAY_L_CURVE
%         figure; 
    subplot(1,2,1); hold on;
    plot(consistency, regularization, 'marker', '.', 'Color',line_color);
    plot(consistency(validIndices==0), regularization(validIndices==0), ...
    'marker', 'x', 'Color','r');
    xlabel('Consistency residual');
    ylabel('Regulariser cost');
    grid on
    box on
    set(gca,'FontSize',9)
end
% Discarding the invalid nodes:
consistency     = consistency(validIndices);
regularization  = regularization(validIndices);
tauList         = tauList(validIndices);        

% cubic spline differentiation to find Kappa (largest curvature) 
eta = log(regularization);
rho = log(consistency);

M = [0 3 0 0;0 0 2 0;0 0 0 1;0 0 0 0];

% spline values for the regularization:
ppE = spline(tauList, eta);
% ppE = polyfit(tauList, eta, 4);  
% keyboard
ppE.coefs = ppE.coefs*M; eta_del  = ppval(ppE,tauList);  
ppE.coefs = ppE.coefs*M; eta_del2 = ppval(ppE,tauList);

% spline values for the consistency:
ppR = spline(tauList, rho);
ppR.coefs = ppR.coefs*M; rho_del  = ppval(ppR,tauList); 
ppR.coefs = ppR.coefs*M; rho_del2 = ppval(ppR,tauList);

% results in the curvature:
Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2).^1.5;

% addpath ~/matlab/ptb/ptb/_Utils
% Kappa = curvature([consistency' regularization'],'polynom',6)   
% Kappa = curvature([rho' eta'],'polynom',4)   
 
% discarding the side nodes: 
% minK = min(Kappa(:)); Kappa([1:2  end-1:end]) = minK; 
Kappa([1:2  end-1:end]) = 0;
% Kappa([1  end]) = 0;

% from which the maximal curvature is at:

index_opt = find(Kappa == max(Kappa), 1, 'first');
disp(['Optimal lambda, Kappa, consistency, regularisation: ', num2str([1/tauList(index_opt), max(Kappa), consistency(index_opt), regularization(index_opt)])])

% disp(['Optimal lambda, Kappa, consistency, regularisation: ', num2str([tauList(index_opt), max(Kappa), consistency(index_opt), regularization(index_opt)])])

% index_opt = find(Kappa == min(Kappa), 1, 'first');
% Kappa = -Kappa;
% disp(['Optimal lambda, Kappa, consistency, regularisation: ', num2str([tauList(index_opt), min(Kappa), consistency(index_opt), regularization(index_opt)])])

if DISPLAY_L_CURVE
%     subplot(1,2,2), semilogx(tauList, Kappa, 'marker', 'x', 'Color',line_color);
%     subplot(1,2,2); hold on; semilogx(tauList(index_opt), Kappa(index_opt), 'marker', 'o', 'Color',line_color);
    subplot(1,2,2), semilogx(1./tauList, Kappa, 'marker', '.', 'Color',line_color);
    subplot(1,2,2); hold on; h1=semilogx(1/tauList(index_opt),Kappa(index_opt),'marker','o','Color',line_color,'MarkerFaceColor','auto');% get(h1)
    subplot(1,2,2); axis([1/(max(tauList(:))) 1/(min(tauList(:))) min(Kappa(:))*1.1 max(Kappa(:))*1.1]);
%     consistency=log(consistency); regularization=log(regularization);
    subplot(1,2,1); hold on; plot(consistency(index_opt), regularization(index_opt),'marker','o','Color',line_color,'MarkerFaceColor','auto');
    subplot(1,2,1); axis([min(consistency(:)) max(consistency(:)) min(regularization(:))*0.9 max(regularization(:))]);
    subplot(1,2,2); 
    xlabel('\lambda');
    ylabel('Curvature');
    grid on
    box on
    set(gca,'FontSize',9)
end

lambda_opt = tauList(index_opt);

% save  Lcurve_validIndices.mat    validIndices
save('Lcurve_idxOpt.txt',       'index_opt',    '-ascii')
save('Lcurve_lambdaOpt.txt',    'lambda_opt',   '-ascii')
save('Lcurve_lambdaList.txt',   'tauList',      '-ascii')
save('Lcurve_Kappa.txt',        'Kappa',        '-ascii')

%     
%     figure; 
%     loglog(consistency, regularization, 'marker', 'x', 'Color',line_color);
%     
