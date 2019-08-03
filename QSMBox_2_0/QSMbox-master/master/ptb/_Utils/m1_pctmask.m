function [pct_msk] = m1_pctmask(ima,msk,pct)
% Returns pct_msk incl. ima's top pct % values within msk
% Created by Julio Acosta-Cabronero, 29 Jun 2017

pu                          = pct/100;
sort_grads                  = sort(ima(msk>0));
nV                          = round(pu*length(sort_grads));
pct_msk                     = zeros(size(ima));
grad_thresh                 = sort_grads(end-nV+1);
pct_msk(ima>=grad_thresh)   = 1;
pct_msk(ima< grad_thresh)   = 0;
