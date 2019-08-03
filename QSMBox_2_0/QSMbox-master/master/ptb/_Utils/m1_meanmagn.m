function [mS,S0,R2star] = m1_meanmagn(magn,TE)
% Synthesise magnitude image for mean TE
%  mS=S0.*exp(-mean(TE).*R2star (S0 & R2star inferred from magn)
% 
% Created by Julio Acosta-Cabronero, 08/12/2018

magn = log(magn);

sz = size(magn);
magn = reshape(magn,[sz(1)*sz(2)*sz(3) sz(4)]);

mTE = mean(TE);
O   = ones([length(TE) 1]);
TE  = [O TE'];

P = TE\magn';
clear magn

S0 = exp(P(1,:)');
R2star = -P(2,:)';
clear P

mS = S0.*exp(-mTE.*R2star);

% S0 = reshape(S0,[sz(1) sz(2) sz(3)]);
% R2star = reshape(R2star,[sz(1) sz(2) sz(3)]);
mS = reshape(mS,[sz(1) sz(2) sz(3)]);
