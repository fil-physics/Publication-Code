function R1_Error_InterScan_MoCo

%
% Assumed settings
%

% Calibration data
TRc = 6.5e-3;
alpha_c = deg2rad(6);

% VFA data
TR1 = 19.5e-3;
TR2 = 19.5e-3;
alpha_1 = deg2rad(6);
alpha_2 = deg2rad(26);

% Empirical knowledge
d12_Empirical = linspace(0.85,1.18,90); % max variance at 7T, d12 = fT1/fT2; max range at 3T: [0.95,1.09]
D12_Empirical = linspace(0.84,1.18,95); % 3T example - assume no B1+ effect, i.e. d = 1


% Simulation settings
nT1 = 60; nfT1 = 80; nfT2 = 100;
T1_range = linspace(0.7,2,nT1);
fT1_range = linspace(0.5,1.5,nfT1);
fT2_range = linspace(0.5,1.5,nfT2);

% [fT1,fT2,R1] = ndgrid(fT1_range,fT2_range,1./T1_range);
[fT1,fT2,R1] = ndgrid(fT1_range,fT2_range,1./T1_range);
alphaE = acos(exp(-TRc.*R1));


%
% dKdD as a function of fT1, fT2 and R1 (via alphaE):
%
dKdD = (fT1.*fT2.^2.*alpha_c^2 + fT1.*alphaE.^2) ...
        ./(fT2.*fT1.^2.*alpha_c.^2 + fT2.*alphaE.^2);

%
% Display
%
% data stored as: fT1 x fT2 x R1 but imagesc maps first axis to vertical 
% i.e. fT1 -> y-axis, fT2 -> x-axis so X and Y passed to imagesc and
% contour are reversed in order...
figure(2); clf; set(gcf, 'Color', [1 1 1]); fontSize = 16;
R1_index = 23;
imagesc(fT2_range,fT1_range,dKdD(:,:,R1_index)); set(gca,'YDir','normal')
hold on; grid on; grid minor
xlabel('f_T_2'); ylabel('f_T_1')
set(gca, 'FontSize',fontSize)

% Conditions insensitive to transmit field effects:
cMap = hot; cVal1 = 250; cVal2 = 60;
fT1_points = linspace(fT1_range(1),fT1_range(end),nfT2);
plot(fT2_range, alphaE(1,:,R1_index).^2./alpha_c.^2./fT1_points, '--', 'Color',cMap(cVal1,:),'LineWidth',2);
plot(alphaE(1,:,R1_index).^2./alpha_c.^2./fT2(1,:,R1_index), fT1_points, '--','Color',cMap(cVal1,:),'LineWidth',2)
plot(fT2_range, fT2_range, '--','Color',cMap(cVal2,:),'LineWidth',2) % fT1 = fT2
axis([fT1_range(1) fT1_range(end) fT2_range(1) fT2_range(end)]) 
legend('f_T_j=\alpha_E^2\alpha_c^{-2}f_T_i^{-1}','f_T_1 = f_T_2','Location','North','box','off')

[c,h] = contour(fT2(:,:,R1_index),fT1(:,:,R1_index),dKdD(:,:,R1_index),0.91:0.06:1.1,'ShowText','on');
h.LevelList = round(h.LevelList,2); h.LineColor = 'k';
clabel(c,h,'FontSize',fontSize); 
h.DisplayName = 'Iso ^{\partial\kappa}/_{\partial\Delta_{1,2}}'; hold off

title(['^{\partial\kappa}/_{\partial\Delta_{1,2}} for R_1 = ' num2str(round(R1(1,1,R1_index)*100)/100) 's^{-1}'])


%**************************************************************************
%
% Compute error in R1 and permute through correction options
%

% d = fT1/fT2 => fT1 = d.*fT2
[d12,D12,R1,fT2] = ndgrid(d12_Empirical,D12_Empirical,1./T1_range,fT2_range);

% Ernst angle
alphaE = acos(exp(-TRc.*R1));

% dKdD as a function of d, fT2 and R1 (via alphaE):
dKdD = ((d12.*fT2).*fT2.^2.*alpha_c^2 + (d12.*fT2).*alphaE.^2) ...
        ./(fT2.*(d12.*fT2).^2.*alpha_c.^2 + fT2.*alphaE.^2);


% Signals:
SPGR = @(R1, TR, alpha) sin(alpha).*(1 - exp(-TR.*R1))./(1 - cos(alpha).*exp(-TR.*R1));

% R1 estimate
smallFA = true;
if smallFA
    R1_est = @(I1, TR1, a1, I2, TR2, a2) 0.5.*(I1.*a1./TR1 - I2.*a2./TR2)./(I2./a2 - I1./a1);
else
    R1_est = @(I1, TR1, a1, I2, TR2, a2) (I1./TR1.*(sin(a2).*(cos(a1)-1)) ...
                                        +I2./TR2.*sin(a1).*(1-cos(a2))) ./ ...
                            ((I1.*cos(a1).*sin(a2)-I2.*cos(a2).*sin(a1)));
end

% VFA data:
% s1 = D12_Empirical(end);s2 = 1; % for simplicity s1 = s2*D12; assume s2 = 1
I1 = D12.*SPGR(R1, TR1, alpha_1.*d12.*fT2); % fT1 = d12.*fT2
I2 = SPGR(R1, TR2, alpha_2.*fT2);

% Time-locked calibration data:
x1 = D12.*SPGR(R1, TRc, alpha_c.*d12.*fT2); % s1 == D12 => s2 = 1; fT1 = d12.*fT2
x2 = SPGR(R1, TRc, alpha_c.*fT2);
K12 = x1./x2;


% Reference R1 under small FA approximation (i.e. small FA error)
% => remove position-specific receive sensitivity effects and use position
% specific transmit field effects (infeasible in practice since s1, s2 unknown):
if smallFA
    R1_Ref = R1_est(I1, TR1, d12.*fT2.*alpha_1, D12.*I2, TR2, fT2.*alpha_2);
else
    R1_Ref = R1;
end

% Typical small FA estimate capturing motion (one Tx field, no Rx corr):
R1_Motion = R1_est(I1, TR1, d12.*fT2.*alpha_1, I2, TR2, d12.*fT2.*alpha_2);
errMotion = (R1_Motion - R1_Ref)./R1_Ref.*100;

% Perfect corrrection of receive effect => only transmit field error:
R1_TxErr = R1_est(I1, TR1, d12.*fT2.*alpha_1, D12.*I2, TR2, d12.*fT2.*alpha_2);
errTx = (R1_TxErr - R1_Ref)./R1_Ref.*100;

% Proposed correction for receive field effects, using only one transmit 
% field (fT1 = d*fT2) measurement.  This introduces an additional transmit 
% field dependence via dKdD (K12 = D12_Empirical(end).*dKdD):
R1_RxCorr = R1_est(I1, TR1, d12.*fT2.*alpha_1, K12.*I2, TR2, d12.*fT2.*alpha_2);
errRxCorr = (R1_RxCorr - R1_Ref)./R1_Ref.*100;

% Perfect per-position corrrection of transmit effect => only receive field
% error. Note: R1_TxByPos == R1_RxErr
R1_TxByPos = R1_est(I1, TR1, d12.*fT2.*alpha_1, I2, TR2, fT2.*alpha_2);
errTxCorr = (R1_TxByPos - R1_Ref)./R1_Ref.*100;

% Accounting for position-specific transmit with receive field correction:
R1_fullCorr = R1_est(I1, TR1, d12.*fT2.*alpha_1, K12.*I2, TR2, fT2.*alpha_2);
errFullCorr = (R1_fullCorr - R1_Ref)./R1_Ref.*100;
                
%
% Display R1 errors:
%
% data stored as: d x fT2 x R1 but imagesc maps first axis to vertical 
% i.e. d -> y-axis, fT2 -> x-axis so X and Y passed to imagesc and
% contour are reversed in order...
if smallFA; figure(3); else figure(4); end
clf; set(gcf, 'Color', [1 1 1]); levels = 3;
fontSize = 14;

subplot 221
imagesc(D12_Empirical,d12_Empirical,errMotion(:,:,R1_index,end/2)); set(gca,'YDir','normal')
hold on; grid on; grid minor
xlabel('\Delta_{1,2}'); ylabel('\delta_{1,2}');
set(gca, 'FontSize',fontSize)
[c,h]=contour(D12(:,:,R1_index,end/2),d12(:,:,R1_index,end/2),errMotion(:,:,R1_index,end/1),levels,'ShowText','on');
h.LevelList = round(h.LevelList,2); h.LineColor = 'k';
clabel(c,h,'FontSize',fontSize); h.DisplayName = 'Error (%)'; hold off
title(['R_1 = ' num2str(round(R1(1,1,R1_index)*100)/100) 's^{-1}.'])
title('Inter-Scan Motion Error')

subplot 222
imagesc(D12_Empirical,d12_Empirical,errTxCorr(:,:,R1_index,end/2)); set(gca,'YDir','normal')
hold on; grid on; grid minor
xlabel('f_T_2'); ylabel('\delta_{1,2}');
xlabel('\Delta_{1,2}'); ylabel('\delta_{1,2}');
set(gca, 'FontSize',fontSize)
[c,h]=contour(D12(:,:,R1_index,end/2),d12(:,:,R1_index,end/2),errTxCorr(:,:,R1_index,end/2),levels,'ShowText','on');
h.LevelList = round(h.LevelList,2); h.LineColor = 'k';
clabel(c,h,'FontSize',fontSize); h.DisplayName = 'Error (%)'; hold off
title(['R_1 = ' num2str(round(R1(1,1,R1_index)*100)/100) 's^{-1}. Tx Corr.'])
title({'Per-position B_1^+'})

subplot 223
imagesc(D12_Empirical,d12_Empirical,errRxCorr(:,:,R1_index,end/2)); set(gca,'YDir','normal')
hold on; grid on; grid minor
xlabel('f_T_2'); ylabel('\delta_{1,2}');
xlabel('\Delta_{1,2}'); ylabel('\delta_{1,2}');
set(gca, 'FontSize',fontSize)
[c,h]=contour(D12(:,:,R1_index,end/2),d12(:,:,R1_index,end/2),errRxCorr(:,:,R1_index,end/2),levels,'ShowText','on');
h.LevelList = round(h.LevelList,2); h.LineColor = 'k';
clabel(c,h,'FontSize',fontSize); h.DisplayName = 'Error (%)'; hold off
title(['R_1 = ' num2str(round(R1(1,1,R1_index)*100)/100) 's^{-1}. Tx Corr.'])
title('Proposed Receive Field Correction')

subplot 224
imagesc(D12_Empirical,d12_Empirical,errFullCorr(:,:,R1_index,end/2)); set(gca,'YDir','normal')
hold on; grid on; grid minor
xlabel('f_T_2'); ylabel('\delta_{1,2}');
xlabel('\Delta_{1,2}'); ylabel('\delta_{1,2}');
set(gca, 'FontSize',fontSize)
[c,h]=contour(D12(:,:,R1_index,end/2),d12(:,:,R1_index,end/2),errFullCorr(:,:,R1_index,end/2),levels,'ShowText','on');
h.LevelList = round(h.LevelList,2); h.LineColor = 'k';
clabel(c,h,'FontSize',fontSize); h.DisplayName = 'Error (%)'; hold off
title(['R_1 = ' num2str(round(R1(1,1,R1_index)*100)/100) 's^{-1}. Full Corr.'])
title('Full Correction')

% Summarise key errors via median of (signed) errors
% Note: some fractions explode when Tx and Rx effects balance such that 
% errMotion ~0 but errTx and errTxCorr are appreciable:
fprintf('Motion can lead to absolute errors up to %3.0f%%\n', max(abs(errMotion(:))))
fprintf('Over the 4D space, a median of %2.0f%% is caused by receive field effects\n', ...
    median(errTxCorr(:)./errMotion(:))*100)
fprintf('While a median of %2.0f%% is caused by transmit field effects\n', ...
    median(errTx(:)./errMotion(:))*100)
