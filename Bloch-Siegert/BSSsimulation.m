function Msig=BSSsimulation(T1,T2,TR,TE,alpha,woff_kHz,B1,T,PhiBaseInc,W1max,W2max,Nspin,Nexc,Fermia,Fermit0,Interleaved)
% Simulation of the imapct of the Bloch-Siegert on the phase of a voxel across TR before and after the BS pulse.
% The voxel is a square grid of Nspin*Nspin isochromats
%
% Inputs:
% - T1 : Longitudinal relaxation time [ms]
% - T2 : Transverse relaxation time [ms]
% - TR : Repetition time [ms]
% - TE : Echo time [ms]
% - alpha : Excitation flip angle [deg]
% - woff_kHz : Off-resonance excitation pulse [kHz]
% - B1 :  Amplitude of the BS pulse [muT]
% - T: Fermi pulse duration [ms]
% - PhiBaseInc: RF spoiling increment [deg]
% - W1max : Crusher dephasing [pi rad]
% - W2max: Spoiler dephasing [pi rad]
% - Nspin: length of the square grid of isochromats
% - Nexc : Number of excitations
% - Fermia : Fermi pulse width
% - Fermit0 : Fermi pulse transition
% - Interleaved : 1 for interleaved acquisition , 0 for sequential acquisition
%
% Outputs:
% - Msig: Matrix (Nex,2) of the phase for every TR:
%               - before the BS pulse and the crushers
%               - after the BS pulse and the crushers
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Fixed parameters
% Gyromagnetic ratio of the proton
gam = 42.6; % Hz/uT
% Maximum phase (constraint of the scanner)
MAXPhase=360000*pi/180;  %rad


%% Initialization

%%% Spoiler dephasing
W2Range = [-(Nspin/2):(Nspin/2-1)]/Nspin*(W2max)*pi;

%%% Crusher dephasing
W1Range= [-(Nspin/2):(Nspin/2-1)]/Nspin*(W1max)*pi;

%%% RF spoiler increment in rad
PhiBaseInc_rad = (PhiBaseInc)/180*pi;


%%% Excitation flip angle in rad
alpha_rad = (alpha)/180*pi;

%%% initialise phase of RF spoiling
Phi = 0;
PhiInc = PhiBaseInc_rad;


%%% Initialize magnetization vector (along z)
M = [zeros(2,Nspin,Nspin);ones(1,Nspin,Nspin)];


%%% Matrix of relaxation
[Ate,Bte] = Relaxation(TE,T1,T2);
[Atr,Btr] = Relaxation(TR-TE,T1,T2);

%%% BS pulse
[BSPulseNorm, time_ms, dDeltaT_ms]=FermiPulseNorm(T,Fermia,Fermit0);
time=time_ms*1e-3; %[s];
dDeltaT=dDeltaT_ms*1e-3; %[s];
BSPulse=B1*BSPulseNorm; %[muT]

Msig = zeros(Nexc,2);
for n=1:Nexc %% loop over excitations
    
    disp(['Excitation pulse number: ' num2str(n) '/' num2str(Nexc)])
    
    if Interleaved==1
        off=n; % reverse sign every excitation via (-1)^off
    else
        if n==(floor(Nexc/2)+1) % reset phase and increment of the RF spoiling after half of the total number of pulses
            Phi=0;
            PhiInc=PhiBaseInc_rad;
            
        end
        if n<(floor(Nexc/2)+1)%reverse sign after half of the pulses
            off=1;  % Negative off-resonance frequency via (-1)^off
        else
            off=2;  % Positive off-resonance frequency via (-1)^off
        end
    end
    
    
    %%% BS pulse  parameters
    woff=((-1)^off)*2*pi*woff_kHz*1000; %% rad/s off resonance frequency of the pulse
    
    % excitation and relaxation
    parfor x=1:Nspin
        A = Ate * Excitation(alpha_rad,Phi);
        B = repmat(Bte,[1 Nspin])
        M(:,:,x) = A*M(:,:,x)+B;
    end
    
    % Signal at TE before BS pulse and crushers
    Msig(n,1) = mean( squeeze(reshape(M(1,:,:),1,[])+1i*reshape(M(2,:,:),1,[]))) * exp(-1i*(Phi));
    
    % BS pulse + crushers
    parfor x=1:Nspin
        for y=1:Nspin
            
            %%%% Crushers and RF Spoiling
            M(:,y,x)=Rotz(-Phi)*Rotz(-W1Range(y))*M(:,y,x);
            
            %%%% BS pulse
            for k=1:length(time) % sampling of the BS pulse
                w1=2*pi*gam*BSPulse(k);
                M(:,y,x)=Rotz(-time(k)*woff)*Roty(w1*dDeltaT)*Rotz(time(k)*woff)*M(:,y,x);
            end
            
            %%%% Crusher with opposite polarity and back to rotating
            %%%% frame of reference
            M(:,y,x)=Rotz(W1Range(y))*Rotz(Phi)*M(:,y,x);
            
        end
    end
    
    
    % Signal at TE after BS pulse
    Msig(n,2) = mean( squeeze((reshape(M(1,:,:),1,[])+1i*reshape(M(2,:,:),1,[]))) ) * exp(-1i*Phi);
    
    % spoiler gradient
    parfor x=1:Nspin
        M(:,:,x) = Rotz(W2Range(x))*squeeze(M(:,:,x));
    end
    
    % relaxation over (TR-TE)
    parfor x=1:Nspin
        M(:,:,x)=Atr*M(:,:,x)+repmat(Btr,[1,Nspin]);
    end
    
    % increment RF phase
    Phi = Phi+PhiInc;
    PhiInc = PhiInc+PhiBaseInc_rad;
    
    Phi=mod(Phi,MAXPhase);
    PhiInc=mod(PhiInc,MAXPhase);
    
end