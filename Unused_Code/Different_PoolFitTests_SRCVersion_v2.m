close all; clear all

TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([14 22 30 38 46 54 62 70]); MaxAngle = 70;
T1_S = 1.15; T1_F = 0.4; T2_S = 0.08; T2_F = 0.02; k_FS = 9; k_SB = 5; k_FB = 5; M0_F = 0.25; M0_S = 0.55; M0_B = 0.2; T1_B = 1; k_SF = (k_FS*M0_F)/M0_S;
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

% GT-signals. Change between CS and US (B1 or TRF) here.
SPGR_Data = CS_SPGR_SteadyState(FA_SPGR,TR_SPGR,MaxAngle,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP0_Data = CS_SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP180_Data = CS_SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data = [SSFP0_Data ; SSFP180_Data];

SNR = 20; Sigma = mean(SPGR_Data)/SNR;

% Add different noise to each element.
SPGR_Data_Noisy = zeros(length(SPGR_Data),1);
for mm = 1:length(SPGR_Data)
    SPGR_Data_Noisy(mm) = SPGR_Data(mm) + (normrnd(0,Sigma));
end
SSFP_Data_Noisy = zeros(length(SSFP_Data),1);
for nn = 1:length(SSFP_Data)
    SSFP_Data_Noisy(nn) = SSFP_Data(nn) + (normrnd(0,Sigma));
end
% Normalise signals and concatenate.
SPGR_Data_Norm = SPGR_Data_Noisy./mean(SPGR_Data_Noisy);
SSFP_Data_Norm = SSFP_Data_Noisy./mean(SSFP_Data_Noisy);
Data_Norm = [SPGR_Data_Norm ; SSFP_Data_Norm];

% Calculate parameter predictions from derived expressions - only valid for CSMT.
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S; dMzB = 0;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
T_RF = deg2rad(MaxAngle)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SSFP);
kSF_Prediction = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
kFS_Prediction = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1F_Prediction = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1S_Prediction = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
M0F_Prediction = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
M0S_Prediction = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
T1F_Prediction = 1/R1F_Prediction; T1S_Prediction = 1/R1S_Prediction;

Plot_SPGR_P = SPGR_SteadyState_nonDE(FA_SPGR,TR_SPGR,'T1_S',T1S_Prediction,'T1_F',T1F_Prediction,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP0_P = SSFP_SteadyState_nonDE(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP180_P = SSFP_SteadyState_nonDE(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP_P = [Plot_SSFP0_P ; Plot_SSFP180_P];

Trials = 5000; Iterations = 7; N = 50; Runs = 1;

% NON-DE FIT.
% [Sol,~,~] = SRC_Sim_NDE_v2(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Norm);
% T1S_final = Sol(1); T1F_final = Sol(2); T2S_final = Sol(3); T2F_final = Sol(4); M0F_final = Sol(5); M0S_final = Sol(6); kFS_final = Sol(7);  kSF_final = Sol(8);
% 
% Plot_SPGR_F = SPGR_SteadyState_nonDE(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SSFP0_F = SSFP_SteadyState_nonDE(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SSFP180_F = SSFP_SteadyState_nonDE(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SSFP_F = [Plot_SSFP0_F ; Plot_SSFP180_F];

% DE FIT.
[Sol,~,~] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Norm);
T1S_final = Sol(1); T1F_final = Sol(2); T2S_final = Sol(3); T2F_final = Sol(4); M0F_final = Sol(5); kFS_final = Sol(6);

Plot_SPGR_F = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'k_FS',kFS_final);
Plot_SSFP0_F = SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'k_FS',kFS_final);
Plot_SSFP180_F = SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'k_FS',kFS_final);
Plot_SSFP_F = [Plot_SSFP0_F ; Plot_SSFP180_F];

figure(1)
subplot(3,1,1)
plot(rad2deg(FA_SPGR), Plot_SPGR_F./mean(Plot_SPGR_F),'LineWidth',2,'LineStyle','--')
hold on
plot(rad2deg(FA_SPGR), SPGR_Data_Norm, 'ko','LineWidth',2)
plot(rad2deg(FA_SPGR), Plot_SPGR_P./mean(Plot_SPGR_P),'LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('SPGR Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12); grid on; grid minor;
subplot(3,1,2)
plot(rad2deg(FA_SSFP0), Plot_SSFP0_F./mean(Plot_SSFP_F),'LineWidth',2,'LineStyle','--')
hold on
plot(rad2deg(FA_SSFP0), SSFP_Data_Norm(1:length(FA_SSFP0)),'ko','LineWidth',2)
plot(rad2deg(FA_SSFP0), Plot_SSFP0_P./mean(Plot_SSFP_P),'LineWidth',2);
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{0} Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12); grid on; grid minor;
ll = legend('Fitted Signal', 'Three-Pool Data','Prediction'); ll.FontSize = 16; legend('boxoff');
subplot(3,1,3)
plot(rad2deg(FA_SSFP180), Plot_SSFP180_F./mean(Plot_SSFP_F),'LineWidth',2,'LineStyle','--')
hold on
plot(rad2deg(FA_SSFP180), SSFP_Data_Norm(length(FA_SSFP0)+1:end),'ko','LineWidth',2)
plot(rad2deg(FA_SSFP180), Plot_SSFP180_P./mean(Plot_SSFP_P),'LineWidth',2);
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{180} Signal','FontSize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12); grid on; grid minor;
