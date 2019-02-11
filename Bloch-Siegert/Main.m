%%  Main.m
%
%   Numerical simulations to evaluate the impact of the Bloch-Siegert pulse on
%	the phase of two acqusiitions with opposite off-resonance frequencies. 
%   The two acquisitions can be obtained in sequential or interleaved order. 
% 
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%%***************************************************************************************%%



%% Define parameters
% Tissue parameters
T1 = 550; %[ms]  Longitudinal relaxation time
T2 = 70;  %[ms]  Transverse relaxation time

% Sequence parameters
TR = 35;  %[ms]  Repetition time
TE = 2;   %[ms]  Echo time
alpha = 20; %[deg]  Excitation flip angle
PhiBaseInc = 137; %[deg] Rf Spoiling increment
W1max= 8; % [pi rad] Crusher gradient dephasing moment
W2max= 4; % [pi rad] Spoiler gradient dephasing moment
Interleaved=1; % 0 for sequential acquisitions , 1 for interleaved acquisition

% Bloch-Siegert pulse parameters
woff_kHz = 2; % [kHz] Off resonance frequency of the BS pulse
B1 = 8; % [uT] BS pulse amplitude 
T = 2; % [ms] BS pulse duration
Fermia = 0.16; % [ms] Fermi pulse "transition" parameter
Fermit0 = 3; %[ms] Fermi pulse "width" parameter

% Simulation parameters
Nspin = 100; % Number of isochromat along each dimension of the 2D grid
Nexc = 100; % Number of excitation pulses simulated

%% Run simulation
Msig=BSSsimulation(T1,T2,TR,TE,alpha,woff_kHz,B1,T,PhiBaseInc,W1max,W2max,Nspin,Nexc,Fermia,Fermit0,Interleaved);


%% Plot results
if Interleaved
    PhiBeforeNeg=angle(Msig(1:2:end,1)); % phase before the BS pulse of the acquisition with negative off-resonance frequency
    PhiAfterNeg=angle(Msig(1:2:end,2)); % phase after the BS pulse of the acquisition with negative off-resonance frequency
    PhiBeforePos=angle(Msig(2:2:end,1)); % phase before the BS pulse of the acquisition with positive off-resonance frequency
    PhiAfterPos=angle(Msig(2:2:end,2));% phase after the BS pulse of the acquisition with positive off-resonance frequency
else
    PhiBeforeNeg=angle(Msig(1:floor(Nexc/2),1));
    PhiAfterNeg=angle(Msig(1:floor(Nexc/2),2));
    PhiBeforePos=angle(Msig((floor(Nexc/2)+1):end,1));
    PhiAfterPos=angle(Msig((floor(Nexc/2)+1):end,2));
end

figure;
subplot(2,1,1)
plot(PhiBeforePos)
hold all
plot(PhiBeforeNeg)
title('Before the BS pulse')
xlabel('Pairs of pulses')
ylabel('Phase [rad]')
legend('-\omega_{off}','+\omega_{off}')

subplot(2,1,2)
plot(PhiAfterPos)
hold all
plot(PhiAfterNeg)
title('After the BS pulse')
xlabel('Pairs of pulses')
ylabel('Phase [rad]')

legend('-\omega_{off}','+\omega_{off}')

