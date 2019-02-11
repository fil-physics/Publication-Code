%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%

function [Fermi, Time, dDeltaT]=FermiPulseNorm(T,Fermia,Fermit0)

% Normalized Fermi pulse as implemented in the MR pulse sequence
% Inputs:
% - T[ms]: duration of the pulse
% - Fermia [ms]:  "transition" parameter
% - Fermit0 [ms]:  "width" parameter
% Ouptuts:
% - Fermi [a.u.]: vector of the normalized Fermi pulse
% - Time [ms]: vector of time stamps
% - dDeltaT [ms]: sampling interval


dDeltaT=0.01; %[ms]  % sampling interval

dTime=-T/2; % left edge of the pulse [ms]

dCorrFactor=T*1e3/8192; % MR Pulse programming requirement

for i=0:(T*1e2)
    Fermi(i+1)=1/(1+exp((abs(dTime)-Fermit0*dCorrFactor)/ (Fermia*dCorrFactor)));
    Time(i+1)=dTime;
    dTime=dTime+dDeltaT;    
end

end
