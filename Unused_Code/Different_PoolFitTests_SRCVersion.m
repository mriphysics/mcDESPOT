close all; clear all

TR_SPGR = 5e-3; TR_SSFP = 5e-3; %TR_SPGR = 5.4e-3; TR_SSFP = 4.4e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([10 13 17 20 23 30 43 60]); FA_SSFP180 = deg2rad([10 13 17 20 23 30 43 60]); MaxAngle = max(rad2deg(FA_SSFP0));
T1_S = 0.965; T1_F = 0.465; T2_S = 0.09; T2_F = 0.012; k_FS = 8; k_SB = 5; k_FB = 5; M0_F = 0.2; M0_S = 0.7; M0_B = 0.1; T1_B = 1; k_SF = (k_FS*M0_F)/M0_S;
PC1 = 0; PC2 = pi;

% GT-signals. Change between CS and US (B1 or TRF) here.
SPGR_Data = CS_SPGR_SteadyState(FA_SPGR,TR_SPGR,MaxAngle,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP0_Data = CS_SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP180_Data = CS_SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP0_Data ; SSFP180_Data];

% Calulcate parameter predictions from derived expressions - only valid for CSMT.
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

Trials = 5000; Iterations = 10; N = 50; Runs = 1;
[T1S_final, T1F_final, M0F_final, M0S_final, kFS_final, kSF_final, T2S_final, T2F_final] = SRC_Sim_NDE(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data);

figure(1)
subplot(3,1,1)
plot(rad2deg(FA_SPGR), (SPGR_SteadyState_nonDE(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SPGR), SPGR_Data, 'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('SPGR Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
str = ['T_{1F}^{app} = ' num2str(round(T1F_final,2)) 's,' ' T_{2F}^{app} = ' num2str(round(T2F_final,2,'significant')) 's,' ' T_{1S}^{app} = ' num2str(round(T1S_final,2)) 's,' ' T_{2S}^{app} = ' num2str(round(T2S_final,2)) 's,' ' M_{0F}^{app} = ' num2str(round(M0F_final,2)) ',' ' M_{0S}^{app} = ' num2str(round(M0S_final,2)) ',' ' k_{FS}^{app} = ' num2str(round(kFS_final)) 's^{-1},' ' k_{SF}^{app} = ' num2str(round(kSF_final)) 's^{-1}']; 
dim = [.695 .12 .205 .06];
an = annotation('textbox',dim,'String',str); an.BackgroundColor = 'w'; an.FaceAlpha = 0.7; an.Color = 'b';

subplot(3,1,2)
plot(rad2deg(FA_SSFP0), (SSFP_SteadyState_nonDE(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP0), SSFP0_Data,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{0} Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
ll = legend('Fitted mcDESPOT Signal', 'Three-Pool Data');
ll.FontSize = 12;
set(ll.BoxFace, 'ColorType','truecoloralpha','ColorData',uint8(255*[.5;.5;.5;.1]));

subplot(3,1,3)
plot(rad2deg(FA_SSFP180), (SSFP_SteadyState_nonDE(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP180), SSFP180_Data,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{180} Signal','FontSize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)

% Test goodness of fit.
% mcDESPOT_Result = [SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final) ; SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final) ; SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)];
% mcDESPOT_Fit = goodnessOfFit(mcDESPOT_Result,Data,'NRMSE');
% % CSMT.
% % Diff_T1F = ((T1F_Prediction - T1F_final)/T1F_Prediction) * 100; Diff_T1S = ((T1S_Prediction - T1S_final)/T1S_Prediction) * 100;
% % Diff_T2F = ((T2_F - T2F_final)/T2_F) * 100; Diff_T2S = ((T2_S - T2S_final)/T2_S) * 100;
% % Diff_kFS = ((kFS_Prediction - kFS_final)/kFS_Prediction) * 100; Diff_kSF = ((kSF_Prediction - kSF_final)/kSF_Prediction) * 100;
% % Diff_M0F = ((M0F_Prediction - M0F_final)/M0F_Prediction) * 100; Diff_M0S = ((M0S_Prediction - M0S_final)/M0S_Prediction) * 100;
% % US
% Diff_T1F = ((T1_F - T1F_final)/T1_F) * 100; Diff_T1S = ((T1_S - T1S_final)/T1_S) * 100;
% Diff_T2F = ((T2_F - T2F_final)/T2_F) * 100; Diff_T2S = ((T2_S - T2S_final)/T2_S) * 100;
% Diff_kFS = ((k_FS - kFS_final)/k_FS) * 100; Diff_kSF = ((k_SF - kSF_final)/k_SF) * 100;
% Diff_M0F = ((M0_F - M0F_final)/M0_F) * 100; Diff_M0S = ((M0_S - M0S_final)/M0_S) * 100;
