%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%


function [Arel,Brel]=Relaxation(T,T1,T2)
%
%	Transformation matrix simulating decay
%	over a time interval T [ms], given relaxation times T1[ms] and T2[ms]

E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Arel = [E2 0 0;0 E2 0;0 0 E1];
Brel = [0 0 1-E1]';