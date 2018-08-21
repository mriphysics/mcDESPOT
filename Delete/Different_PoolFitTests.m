%% Fitting a single-pool and mcDESPOT model to 3-pool US and CSMT data. Make sure 2-pool signal functions fit for kSF and M0S also.

clear all; close all

Runs = 10; options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8);

%TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]); MaxAngle = 65;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]); MaxAngle = 70;

T1_S = 1.15; T1_F = 0.4; T2_S = 0.08; T2_F = 0.02; k_FS = 9; k_SB = 5; k_FB = 5; M0_F = 0.25; M0_S = 0.55; M0_B = 0.2; T1_B = 1; k_SF = (k_FS*M0_F)/M0_S;
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

% GT-signals.
SPGR_Data_CS = CS_SPGR_SteadyState(FA_SPGR, TR_SPGR, MaxAngle,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP0_Data_CS = CS_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP180_Data_CS = CS_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_CS = [SSFP0_Data_CS ; SSFP180_Data_CS];
%Data_CS = [SPGR_Data_CS./mean(SPGR_Data_CS) ; SSFP0_Data_CS./mean(SSFP_Data_CS) ; SSFP180_Data_CS./mean(SSFP_Data_CS)];
%FA_Vector = rad2deg([FA_SPGR , FA_SSFP0 , FA_SSFP180]');

SPGR_Data_B1 = B1_SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP0_Data_B1 = B1_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP180_Data_B1 = B1_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_B1 = [SSFP0_Data_B1 ; SSFP180_Data_B1];
%Data_B1 = [SPGR_Data_B1./mean(SPGR_Data_B1) ; SSFP0_Data_B1./mean(SSFP_Data_B1) ; SSFP180_Data_B1./mean(SSFP_Data_B1)];

SPGR_Data_TRF = TRF_SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP0_Data_TRF = TRF_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP180_Data_TRF = TRF_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_TRF = [SSFP0_Data_TRF ; SSFP180_Data_TRF];
%Data_TRF = [SPGR_Data_TRF./mean(SPGR_Data_TRF) ; SSFP0_Data_TRF./mean(SSFP_Data_TRF) ; SSFP180_Data_TRF./mean(SSFP_Data_TRF)];

%% Predictions.

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

%% mcDESPOT fit.

% T1S_LB = 0.25; T1S_UB = 1.5; T2S_LB = 0.04; T2S_UB = 0.15;
% T1F_LB = 0.1; T1F_UB = 0.8; T2F_LB = 0.01; T2F_UB = 0.03;
% M0F_LB = 0.001; M0F_UB = 0.35; M0S_LB = 0.001; M0S_UB = 0.8;
% kFS_LB = (1/0.6); kFS_UB = 40; kSF_LB = (1/0.6); kSF_UB = 20;
% Delta_LB = -pi; Delta_UB = pi; % Wood.
% 
% T1S_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); T2F_app = zeros(Runs,1); M0F_app = zeros(Runs,1); M0S_app = zeros(Runs,1); kFS_app = zeros(Runs,1); kSF_app = zeros(Runs,1); Delta_app = zeros(Runs,1);
% T1S_rand = zeros(Runs,1); T2S_rand = zeros(Runs,1);
% T1F_rand = zeros(Runs,1); T2F_rand = zeros(Runs,1);
% M0F_rand = zeros(Runs,1); kFS_rand = zeros(Runs,1);
% M0S_rand = zeros(Runs,1); kSF_rand = zeros(Runs,1);
% Delta_rand = zeros(Runs,1);
% 
% Fval = zeros(Runs,1); Exitflag = zeros(Runs,1);
% 
% tic
% for ii = 1:Runs
%     
%     disp(ii)
%     
%     % Randomly specify parameter values for fitting starting point.
%     T1S_rand(ii) = T1S_LB + (T1S_UB - T1S_LB) .* rand(1,1); T2S_rand(ii) = T2S_LB + (T2S_UB - T2S_LB) .* rand(1,1);
%     T1F_rand(ii) = T1F_LB + (T1F_UB - T1F_LB) .* rand(1,1); T2F_rand(ii) = T2F_LB + (T2F_UB - T2F_LB) .* rand(1,1);
%     M0F_rand(ii) = M0F_LB + (M0F_UB - M0F_LB) .* rand(1,1); kFS_rand(ii) = kFS_LB + (kFS_UB - kFS_LB) .* rand(1,1);
%     M0S_rand(ii) = M0S_LB + (M0S_UB - M0S_LB) .* rand(1,1); kSF_rand(ii) = kSF_LB + (kSF_UB - kSF_LB) .* rand(1,1);
%     Delta_rand(ii) = Delta_LB + (Delta_UB - Delta_LB) .* rand(1,1);
%     x0 = [T1S_rand(ii) T1F_rand(ii) M0F_rand(ii) M0S_rand(ii) kFS_rand(ii) kSF_rand(ii) T2S_rand(ii) T2F_rand(ii) Delta_rand(ii)];
%     
%     Sig_SPGR = @(x)(SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig_SSFP0 = @(x)(SSFP_SteadyState(FA_SSFP0, TR_SSFP, (0+x(9)),'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig_SSFP180 = @(x)(SSFP_SteadyState(FA_SSFP180, TR_SSFP, (pi+x(9)),'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig_SSFP = @(x)[Sig_SSFP0(x) ; Sig_SSFP180(x)];
%     Sig = @(x)([Sig_SPGR(x)./mean(Sig_SPGR(x)) ; Sig_SSFP(x)./mean(Sig_SSFP(x))]);
%     CF = @(x)sum((Sig(x) - Data).^2);
%     [Sol,Fval(ii),Exitflag(ii),~,~,~,~] = fmincon(CF,x0,[],[],[],[],[T1S_LB T1F_LB M0F_LB M0S_LB kFS_LB kSF_LB T2S_LB T2F_LB Delta_LB],[T1S_UB T1F_UB M0F_UB M0S_UB kFS_UB kSF_UB T2S_UB T2F_UB Delta_UB], [], options);
%     
%     T1S_app(ii) = Sol(1); T1F_app(ii) = Sol(2); M0F_app(ii) = Sol(3); M0S_app(ii) = Sol(4); kFS_app(ii) = Sol(5); kSF_app(ii) = Sol(6); T2S_app(ii) = Sol(7); T2F_app(ii) = Sol(8); Delta_app(ii) = Sol(9);
%     
% end
% toc
% [~,Index] = min(Fval);
% T1S_final = T1S_app(Index); T1F_final = T1F_app(Index); M0F_final = M0F_app(Index); M0S_final = M0S_app(Index); kFS_final = kFS_app(Index); kSF_final = kSF_app(Index); T2S_final = T2S_app(Index); T2F_final = T2F_app(Index); Delta_final = Delta_app(Index);

%% Signal plots. Change title and subplot indices for CS or US.

% Plot_SPGR = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SPGR_P = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_Prediction,'T1_F',T1F_Prediction,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
% Plot_SSFP0 = SSFP_SteadyState(FA_SSFP0,TR_SSFP,Delta_final,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SSFP0_P = SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
% Plot_SSFP180 = SSFP_SteadyState(FA_SSFP180,TR_SSFP,pi+Delta_final,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final);
% Plot_SSFP180_P = SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
% Plot_SSFP = [Plot_SSFP0 ; Plot_SSFP180];
% Plot_SSFP_P = [Plot_SSFP0_P ; Plot_SSFP180_P];
% 
% figure(1)
% subplot(3,2,2)
% plot(rad2deg(FA_SPGR), SPGR_Data./mean(SPGR_Data), 'bo','LineWidth',2); hold on
% plot(rad2deg(FA_SPGR), Plot_SPGR./mean(Plot_SPGR),'k-','LineWidth',2);
% plot(rad2deg(FA_SPGR), Plot_SPGR_P./mean(Plot_SPGR_P),'r--','LineWidth',2);
% xlabel('FA [deg]','FontSize',12); ylabel('SPGR Signal','Fontsize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
% str = ['T_{1F}^{app} = ' num2str(round(T1F_final,2)) 's,' ' T_{2F}^{app} = ' num2str(round(T2F_final,2,'significant')) 's,' ' T_{1S}^{app} = ' num2str(round(T1S_final,2)) 's,' ' T_{2S}^{app} = ' num2str(round(T2S_final,2)) 's,' ' M_{0F}^{app} = ' num2str(round(M0F_final,2)) ',' ' M_{0S}^{app} = ' num2str(round(M0S_final,2)) ',' ' k_{FS}^{app} = ' num2str(round(kFS_final)) 's^{-1},' ' k_{SF}^{app} = ' num2str(round(kSF_final)) 's^{-1}']; 
% dim = [.690 .12 .21 .06];
% an = annotation('textbox',dim,'String',str); an.BackgroundColor = 'w';
% title('CSMT', 'FontSize', 14)
% 
% subplot(3,2,4)
% plot(rad2deg(FA_SSFP0), SSFP0_Data./mean(SSFP0_Data),'bo','LineWidth',2); hold on
% plot(rad2deg(FA_SSFP0), Plot_SSFP0./mean(Plot_SSFP),'k-','LineWidth',2);
% plot(rad2deg(FA_SSFP0), Plot_SSFP0_P./mean(Plot_SSFP_P),'r--','LineWidth',2);
% xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{0} Signal','Fontsize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
% ll = legend('Three-Pool Data','Fitted mcDESPOT Signal','Predicted Signal');
% ll.FontSize = 12;
% 
% subplot(3,2,6)
% plot(rad2deg(FA_SSFP180), SSFP180_Data./mean(SSFP180_Data),'bo','LineWidth',2); hold on
% plot(rad2deg(FA_SSFP180), Plot_SSFP180./mean(Plot_SSFP),'k-','LineWidth',2);
% plot(rad2deg(FA_SSFP180), Plot_SSFP180_P./mean(Plot_SSFP_P),'r--','LineWidth',2);
% xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{180} Signal','FontSize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)

Plot_SPGR_P = SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_Prediction,'T1_F',T1F_Prediction,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP0_P = SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP180_P = SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1S_Prediction,'T2_S',T2_S,'T1_F',T1F_Prediction,'T2_F',T2_F,'M0_F',M0F_Prediction,'M0_S',M0S_Prediction,'k_FS',kFS_Prediction,'k_SF',kSF_Prediction);
Plot_SSFP_P = [Plot_SSFP0_P ; Plot_SSFP180_P];

figure(1);
P1 = plot(rad2deg(FA_SPGR), SPGR_Data_CS./mean(SPGR_Data_CS), 'mo','LineWidth',2); hold on
P2 = plot(rad2deg(FA_SSFP0), SSFP0_Data_CS./mean(SSFP_Data_CS),'mo','LineWidth',2);
P3 = plot(rad2deg(FA_SSFP180), SSFP180_Data_CS./mean(SSFP_Data_CS),'mo','LineWidth',2);

P4 = plot(rad2deg(FA_SPGR), SPGR_Data_B1./mean(SPGR_Data_B1), 'b:','LineWidth',2);
P5 = plot(rad2deg(FA_SSFP0), SSFP0_Data_B1./mean(SSFP_Data_B1),'b:','LineWidth',2);
P6 = plot(rad2deg(FA_SSFP180), SSFP180_Data_B1./mean(SSFP_Data_B1),'b:','LineWidth',2);

P7 = plot(rad2deg(FA_SPGR), SPGR_Data_TRF./mean(SPGR_Data_TRF), 'r-.','LineWidth',2);
P8 = plot(rad2deg(FA_SSFP0), SSFP0_Data_TRF./mean(SSFP_Data_TRF),'r-.','LineWidth',2);
P9 = plot(rad2deg(FA_SSFP180), SSFP180_Data_TRF./mean(SSFP_Data_TRF),'r-.','LineWidth',2);

P10 = plot(rad2deg(FA_SPGR), Plot_SPGR_P./mean(Plot_SPGR_P),'k-','LineWidth',2); 
P11 = plot(rad2deg(FA_SSFP0), Plot_SSFP0_P./mean(Plot_SSFP_P),'k-','LineWidth',2);
P12 = plot(rad2deg(FA_SSFP180), Plot_SSFP180_P./mean(Plot_SSFP_P),'k-','LineWidth',2);

ll = legend([P1 P4 P7 P10],{'CSMT','US (B_{1}-Scaling)','US (T_{RF}-Scaling)','Predicted'});
ll.FontSize = 18; legend boxoff
xlabel('FA (deg)', 'FontSize', 18); ylabel('Normalised Signal (a.u.)', 'FontSize', 18);
get(gca, 'XTick');
set(gca, 'FontSize', 16)
get(gca, 'YTick');
set(gca, 'FontSize', 16)
grid on; grid minor;

%% Outtakes.

% Test goodness of fit.
% SP_Result = [SPGR_SP_SteadyState(FA_SPGR,TR_SPGR,'T1',T1_final,'M0',M0_final) ; SSFP_SP_SteadyState(FA_SSFP,TR_SSFP,PC1,'T1',T1_final,'T2',T2_final,'M0',M0_final) ; SSFP_SP_SteadyState(FA_SSFP,TR_SSFP,PC2,'T1',T1_final,'T2',T2_final,'M0',M0_final)];
% mcDESPOT_Result = [SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final) ; SSFP_SteadyState(FA_SSFP,TR_SSFP,PC1,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final) ; SSFP_SteadyState(FA_SSFP,TR_SSFP,PC2,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)];
% SP_Fit  = goodnessOfFit(SP_Result, Data, 'NRMSE');
% mcDESPOT_Fit = goodnessOfFit(mcDESPOT_Result,Data,'NRMSE');
% CSMT.
% Diff_T1F = ((T1F_Prediction - T1F_final)/T1F_Prediction) * 100; Diff_T1S = ((T1S_Prediction - T1S_final)/T1S_Prediction) * 100;
% Diff_T2F = ((T2_F - T2F_final)/T2_F) * 100; Diff_T2S = ((T2_S - T2S_final)/T2_S) * 100;
% Diff_kFS = ((kFS_Prediction - kFS_final)/kFS_Prediction) * 100; Diff_kSF = ((kSF_Prediction - kSF_final)/kSF_Prediction) * 100;
% Diff_M0F = ((M0F_Prediction - M0F_final)/M0F_Prediction) * 100; Diff_M0S = ((M0S_Prediction - M0S_final)/M0S_Prediction) * 100;
% US
% Diff_T1F = ((T1_F - T1F_final)/T1_F) * 100; Diff_T1S = ((T1_S - T1S_final)/T1_S) * 100;
% Diff_T2F = ((T2_F - T2F_final)/T2_F) * 100; Diff_T2S = ((T2_S - T2S_final)/T2_S) * 100;
% Diff_kFS = ((k_FS - kFS_final)/k_FS) * 100; Diff_kSF = ((k_SF - kSF_final)/k_SF) * 100;
% Diff_M0F = ((M0_F - M0F_final)/M0_F) * 100; Diff_M0S = ((M0_S - M0S_final)/M0_S) * 100;

% Fitting a mcDESPOT model to 2-pool MT data.

% Runs = 10;
% options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8);
% 
% % Define sequence and ground-truth tissue parameters.
% % TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
% % FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
% %T1_W = 1.15; T2_W = 0.08; M0_B = 0.25; k_WB = 9; T1_B = 0.4;
% 
% TR_SPGR = 5.4e-3; TR_SSFP = 4.4e-3;
% FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([10 13 17 20 23 30 43 60]);
% T1_W = 1.15; T2_W = 0.08; M0_B = 0.15; k_WB = 9; T1_B = 1;
% 
% T1S_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); T2F_app = zeros(Runs,1); M0F_app = zeros(Runs,1); M0S_app = zeros(Runs,1); kFS_app = zeros(Runs,1); kSF_app = zeros(Runs,1);
% T1S_rand = zeros(Runs,1); T2S_rand = zeros(Runs,1);
% T1F_rand = zeros(Runs,1); T2F_rand = zeros(Runs,1);
% M0F_rand = zeros(Runs,1); kFS_rand = zeros(Runs,1);
% M0S_rand = zeros(Runs,1); kSF_rand = zeros(Runs,1);
% 
% Fval = zeros(Runs,1); Exitflag = zeros(Runs,1); Grad = zeros(8, Runs); Hessian = zeros(8, 8, Runs);
% 
% % GT-signals. Change between CS and US in functions.
% SPGR_Data = TwoPoolMT_SPGR(FA_SPGR,TR_SPGR,'T1_W',T1_W,'T1_B',T1_B,'M0_B',M0_B,'k_WB',k_WB);
% SSFP_Data = TwoPoolMT_SSFP(FA_SSFP,TR_SSFP,'T1_W',T1_W,'T1_B',T1_B,'T2_W',T2_W,'M0_B',M0_B,'k_WB',k_WB);
% SSFP_Data_180 = TwoPoolMT_SSFP_180(FA_SSFP,TR_SSFP,'T1_W',T1_W,'T1_B',T1_B,'T2_W',T2_W,'M0_B',M0_B,'k_WB',k_WB);
% Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];
% 
% T1S_LB = 0.2; T1S_UB = 1.25; T2S_LB = 0.04; T2S_UB = 0.15;
% T1F_LB = 0.1; T1F_UB = 0.5; T2F_LB = 0.01; T2F_UB = 0.03;
% M0F_LB = 0.05; M0F_UB = 0.3; M0S_LB = 0.05; M0S_UB = 0.8;
% kFS_LB = 1; kFS_UB = 20; kSF_LB = 1; kSF_UB = 10; % Original bounds.
% 
% % T1S_LB = 0; T1S_UB = 3; T2S_LB = 0; T2S_UB = 0.5;
% % T1F_LB = 0; T1F_UB = 1; T2F_LB = 0; T2F_UB = 0.3;
% % M0F_LB = 0; M0F_UB = 0.5; M0S_LB = 0; M0S_UB = 0.8;
% % kFS_LB = 0; kFS_UB = 40; kSF_LB = 0; kSF_UB = 20;
% 
% tic
% for ii = 1:Runs
%     
%     disp(ii)
%     
%     % Randomly specify parameter values for fitting starting point.
%     T1S_rand(ii) = T1S_LB + (T1S_UB - T1S_LB) .* rand(1,1); T2S_rand(ii) = T2S_LB + (T2S_UB - T2S_LB) .* rand(1,1);
%     T1F_rand(ii) = T1F_LB + (T1F_UB - T1F_LB) .* rand(1,1); T2F_rand(ii) = T2F_LB + (T2F_UB - T2F_LB) .* rand(1,1);
%     M0F_rand(ii) = M0F_LB + (M0F_UB - M0F_LB) .* rand(1,1); kFS_rand(ii) = kFS_LB + (kFS_UB - kFS_LB) .* rand(1,1);
%     M0S_rand(ii) = M0S_LB + (M0S_UB - M0S_LB) .* rand(1,1); kSF_rand(ii) = kSF_LB + (kSF_UB - kSF_LB) .* rand(1,1);
%     x0 = [T1S_rand(ii) T1F_rand(ii) M0F_rand(ii) M0S_rand(ii) kFS_rand(ii) kSF_rand(ii) T2S_rand(ii) T2F_rand(ii)];
%     
%     Sig_SPGR = @(x)(SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig_SSFP = @(x)(SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig_SSFP180 = @(x)(SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
%     Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
%     CF = @(x)norm(Sig(x) - Data).^2;
%     [Sol,Fval(ii),Exitflag(ii),Output,Lambda,Grad(:,ii),Hessian(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1S_LB T1F_LB M0F_LB M0S_LB kFS_LB kSF_LB T2S_LB T2F_LB],[T1S_UB T1F_UB M0F_UB M0S_UB kFS_UB kSF_UB T2S_UB T2F_UB], [], options);
%     
%     % Remove non-convergent solutions.
%     if Exitflag(ii)  == 1
%         
%         T1S_app(ii) = Sol(1); T1F_app(ii) = Sol(2); M0F_app(ii) = Sol(3); M0S_app(ii) = Sol(4); kFS_app(ii) = Sol(5); kSF_app(ii) = Sol(6); T2S_app(ii) = Sol(7); T2F_app(ii) = Sol(8);
%         
%     end
%     
% end
% toc
% 
% [~,Index] = min(Fval);
% 
% T1S_final = T1S_app(Index); T1F_final = T1F_app(Index); M0F_final = M0F_app(Index); M0S_final = M0S_app(Index); kFS_final = kFS_app(Index); kSF_final = kSF_app(Index); T2S_final = T2S_app(Index); T2F_final = T2F_app(Index);
% 
% % Signal plots. Change title and subplot indices for CS or US.
% 
% figure(1)
% subplot(3,3,1)
% plot(rad2deg(FA_SPGR), (SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
% hold on
% plot(rad2deg(FA_SPGR), SPGR_Data, 'ro','LineWidth',2)
% xlabel('FA [deg]','FontSize',12); ylabel('SPGR Signal','Fontsize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
% title('B_{1}-Scaling','FontSize',14);
% str = ['T_{1F}^{app} = ' num2str(round(T1F_final,2)) 's,' ' T_{2F}^{app} = ' num2str(round(T2F_final,2,'significant')) 's,' ' T_{1S}^{app} = ' num2str(round(T1S_final,2)) 's,' ' T_{2S}^{app} = ' num2str(round(T2S_final,2)) 's,' ' M_{0F}^{app} = ' num2str(round(M0F_final,2)) ',' ' M_{0S}^{app} = ' num2str(round(M0S_final,2)) ',' ' k_{FS}^{app} = ' num2str(round(kFS_final)) 's^{-1},' ' k_{SF}^{app} = ' num2str(round(kSF_final)) 's^{-1}']; 
% dim = [.14 .42 .2 .075];
% %dim = [.42 .42 .2 .075];
% %dim = [.7 .42 .2 .075];
% an = annotation('textbox',dim,'String',str); an.BackgroundColor = 'w'; an.FaceAlpha = 0.8;
% 
% subplot(3,3,7)
% plot(rad2deg(FA_SSFP), (SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
% hold on
% plot(rad2deg(FA_SSFP), SSFP_Data,'ro','LineWidth',2)
% xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{0} Signal','Fontsize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
% ll = legend('Fitted mcDESPOT Signal', 'Two-Pool MT Data');
% ll.FontSize = 12;
% set(ll.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.5;.5;.5;.1]));
% 
% subplot(3,3,4)
% plot(rad2deg(FA_SSFP), (SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
% hold on
% plot(rad2deg(FA_SSFP), SSFP_Data_180,'ro','LineWidth',2)
% xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{180} Signal','FontSize',12)
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)


% % Single-pool fit.
% %T1_LB = 0.1; T1_UB = 1.25; T2_LB = 0.01; T2_UB = 0.15; M0_LB = 0.05; M0_UB = 1.5;
% T1_LB = 0.1; T1_UB = 1.5; T2_LB = 0.01; T2_UB = 0.15; M0_LB = 0.05; M0_UB = 1;
% T1_app = zeros(Runs,1); T2_app = zeros(Runs,1); M0_app = zeros(Runs,1);
% T1_rand = zeros(Runs,1); T2_rand = zeros(Runs,1); M0_rand = zeros(Runs,1);
% 
% Fval_SP = zeros(Runs,1); Exitflag_SP = zeros(Runs,1); Grad_SP = zeros(3, Runs); Hessian_SP = zeros(3, 3, Runs);
% 
% tic
% for ii = 1:Runs
%     
%     disp(ii)
%     
%     % Randomly specify parameter values for fitting starting point.
%     T1_rand(ii) = T1_LB + (T1_UB - T1_LB) .* rand(1,1); T2_rand(ii) = T2_LB + (T2_UB - T2_LB) .* rand(1,1);
%     M0_rand(ii) = M0_LB + (M0_UB - M0_LB) .* rand(1,1);
%     x0 = [T1_rand(ii) M0_rand(ii) T2_rand(ii)];
%     
%     Sig_SPGR = @(x)(SPGR_SP_SteadyState(FA_SPGR,TR_SPGR,'T1',x(1),'M0',x(2)));
%     Sig_SSFP = @(x)(SSFP_SP_SteadyState(FA_SSFP,TR_SSFP,0,'T1',x(1),'T2',x(3),'M0',x(2)));
%     Sig_SSFP180 = @(x)(SSFP_SP_SteadyState(FA_SSFP,TR_SSFP,pi,'T1',x(1),'T2',x(3),'M0',x(2)));
%     Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
%     CF = @(x)sum((Sig(x) - Data).^2);
%     [Sol_SP,Fval_SP(ii),Exitflag_SP(ii),Output,Lambda,Grad_SP(:,ii),Hessian_SP(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1_LB M0_LB T2_LB],[T1_UB M0_UB T2_UB], [], options);
%     
%     % Remove non-convergent solutions.
%     if Exitflag_SP(ii) == 1
%         
%         T1_app(ii) = Sol_SP(1); M0_app(ii) = Sol_SP(2); T2_app(ii) = Sol_SP(3);
%     
%     end
%     
% end
% toc
% [~,Index] = min(Fval_SP);
% T1_final = T1_app(Index); M0_final = M0_app(Index); T2_final = T2_app(Index);