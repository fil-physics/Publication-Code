function lambda_opt = m1_Lcurve(lambda_list_fname,scales)
% DESCRIPTION
%  Estimate optimal regularisation parameter using L-curve analysis
% 
% SINTAX
%  lambda_opt = m1_Lcurve(lambda_list,scales)
% 
% INPUTS
%  lambda_list_fname    Text file containing a list of processed lambdas
% 
%  scales               Vector of scale(s) to be analysed, e.g. 1:4. For
%                       nMEDI, nMEDI with SMV sharpening or if only 
%                       interested in one MSDI scale, 'scales' can be a 
%                       single integer
% 
% Created by Julio Acosta-Cabronero

close all
figure

lambda_list = load(lambda_list_fname);
lambda_list = lambda_list([2:2:9 11:end]);
disp(num2str(lambda_list))

% c=['m','b','r','k','g'];
c=['k','k','k','k','k'];

for S=scales
    for x=1:length(lambda_list)
%        fid_cost_ttt = load(['cost_fid_170706_l' num2str(lambda_list(x)) '.txt']);
%        reg_cost_ttt = load(['cost_reg_170706_l' num2str(lambda_list(x)) '.txt']);
         fid_cost_ttt = load(['cost_fid_l' num2str(lambda_list(x)) '.txt']);
         reg_cost_ttt = load(['cost_reg_l' num2str(lambda_list(x)) '.txt']);
        fid_cost_ttt = sum([fid_cost_ttt(1:S)].^1);
        reg_cost_ttt = sum(reg_cost_ttt(1:S));
        idx=find(fid_cost_ttt>0);
        fid_cost(x)=fid_cost_ttt(idx(end));
        reg_cost(x)=reg_cost_ttt(idx(end));
    end

    [idx_opt,lambda_opt_ttt,Kappa] = m1_LcurveOpt(fid_cost,reg_cost,1./lambda_list,c(S));
    lambda_opt(S) = 1/lambda_opt_ttt;
    
%     [idx_opt,lambda_opt_ttt,Kappa] = m1_LcurveOpt(fid_cost,reg_cost,lambda_list,c(S));
%     lambda_opt(S) = lambda_opt_ttt;
    
end