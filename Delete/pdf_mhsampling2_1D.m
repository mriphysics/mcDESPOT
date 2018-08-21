%% Generates heat maps to depict cost-function space. MaxAngle is in degrees and phase-cycling is in radians.

close all; clear all

% Tissue and sequence parameters.
TR_SPGR = 5e-3; TR_SSFP = 5e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([10 13 17 20 23 30 43 60]); FA_SSFP180 = deg2rad([10 13 17 20 23 30 43 60]); MaxAngle = 60;
T1_S = 0.965; T1_F = 0.465; T1_B = 1; T2_S = 0.09; T2_F = 0.012; M0_F = 0.2; k_FS = 8; M0_S = 0.7; M0_B = 0.1; k_SB = 5; k_FB = 5; k_SF = (M0_F*k_FS)/M0_S; 
PC1 = 0; PC2 = pi; Sigma = 1/300; %SNR = 30.

% TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; 
% FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([2 5 10 15 20 30 40 50]); FA_SSFP180 = deg2rad([2 5 10 15 20 30 40 50]); MaxAngle = 50;
% T1_S = 1.15; T1_F = 0.4; T1_B = 1; T2_S = 0.08; T2_F = 0.02; M0_F = 0.25; M0_S = 0.55; M0_B = 0.2; k_FS = 9; k_SB = 5; k_FB = 5; k_SF = (k_FS*M0_F)/M0_S;
% PC1 = 0; PC2 = pi; Sigma = 1/300; %SNR = 30.

% Add noise to parameters - unnecessary?
% Parameters = [T1_S ; T1_F ; T1_B ; T2_S ; T2_F ; M0_F ; k_FS ; M0_S ; M0_B ; k_SB ; k_FB]; Sigma = 1/300; 
% Parameters_Noisy = zeros(length(Parameters),1); 
% for ii = 1:length(Parameters)
%     Parameters_Noisy(ii) = Parameters(ii) + (normrnd(0,Sigma)); % Scale as T2F can go negative?
% end
% T1_S = Parameters_Noisy(1); T1_F = Parameters_Noisy(2); T1_B = Parameters_Noisy(3); 
% T2_S = Parameters_Noisy(4); T2_F = Parameters_Noisy(5);
% M0_F = Parameters_Noisy(6); k_FS = Parameters_Noisy(7);
% M0_S = Parameters_Noisy(8); M0_B = Parameters_Noisy(9); 
% k_SB = Parameters_Noisy(10); k_FB = Parameters_Noisy(11); 

%% Ground-truth signals. 

% CSMT GT.
SPGR_Data = TRF_SPGR_SteadyState(FA_SPGR, TR_SPGR, MaxAngle,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_0 = TRF_SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = TRF_SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2, MaxAngle,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

% Add noise to data.
% Data_Noisy = zeros(length(Data),1); Sigma_Data = 1/300;
% for jj = 1:length(Data)
%     Data_Noisy(jj) = Data(jj) + (normrnd(0,Sigma_Data));
% end

Upper = [1.25 0.5 0.30 0.80 20 10 0.15 0.03]; Lower = [0.20 0.1 0.05 0.05 1 1 0.04 0.01];

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

%% Direct sampling of cost-function in 2D.

Steps = 100; nTrials = 20000;
T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); M0S_Vector = linspace(Lower(4),Upper(4),Steps);
kFS_Vector = linspace(Lower(5),Upper(5),Steps); kSF_Vector = linspace(Lower(6),Upper(6),Steps);
T2S_Vector = linspace(Lower(7),Upper(7),Steps); T2F_Vector = linspace(Lower(8),Upper(8),Steps);

P_T1S = zeros(length(T1S_Vector),1); P_T1F = zeros(length(T1S_Vector),1);
P_T2S = zeros(length(T2S_Vector),1); P_T2F = zeros(length(T2S_Vector),1);
P_kFS = zeros(length(kFS_Vector),1); P_kSF = zeros(length(kFS_Vector),1);
P_M0F = zeros(length(M0F_Vector),1); P_M0S = zeros(length(M0F_Vector),1);

for ii = 1:length(T1S_Vector)
    tic
    disp(['Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    %%% CHANGE INPUT DATA TO DATA_NOISY? %%%
    for nn = 1:nTrials
        
        %if nn == 1
        %
        %T1F_Rand = T1F_Prediction; T1S_Rand = T1S_Prediction;
        %T2F_Rand = T2_F; T2S_Rand = T2_S;
        %kFS_Rand = kFS_Prediction; kSF_Rand = kSF_Prediction;
        %M0F_Rand = M0F_Prediction; M0S_Rand = M0S_Prediction;
        %
        %else
        
        % Draw random parameter values from a uniform distribution.
        T1S_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
        T1F_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
        M0F_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
        M0S_Rand = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
        kFS_Rand = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
        kSF_Rand = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);
        T2S_Rand = (Upper(7) - Lower(7)) .* rand(1,1) + Lower(7);
        T2F_Rand = (Upper(8) - Lower(8)) .* rand(1,1) + Lower(8);
        
        %end
        
        P_T1S(ii,nn) = (logpdf([T1S_Vector(ii) T1F_Rand M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1F(ii,nn) = (logpdf([T1S_Rand T1F_Vector(ii) M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Vector(ii) T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Vector(ii)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Rand kFS_Vector(ii) kSF_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kSF(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Rand kFS_Rand kSF_Vector(ii) T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_M0F(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Vector(ii) M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
        P_M0S(ii,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Vector(ii) kFS_Rand kSF_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
    end
    toc
end

Max_T1S = max(exp(P_T1S),[],2); Max_T1F = max(exp(P_T1F),[],2);
Max_T2S = max(exp(P_T2S),[],2); Max_T2F = max(exp(P_T2F),[],2);
Max_kFS = max(exp(P_kFS),[],2); Max_kSF = max(exp(P_kSF),[],2);
Max_M0F = max(exp(P_M0F),[],2); Max_M0S = max(exp(P_M0S),[],2);

subplot(2,4,1); bar(T1F_Vector, Max_T1F, 'FaceColor', 'k'); xlabel('T_{1F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(T1_F, 'r--'); xlim([Lower(2) Upper(2)]); %ylim([0.4 1]);
subplot(2,4,2); bar(T1S_Vector, Max_T1S, 'FaceColor', 'k'); xlabel('T_{1S}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(T1_S, 'r--'); xlim([Lower(1) Upper(1)]); %ylim([0.4 1]);
subplot(2,4,3); bar(T2F_Vector, Max_T2F, 'FaceColor', 'k'); xlabel('T_{2F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(T2_F, 'r--'); xlim([Lower(8) Upper(8)]); %ylim([0.4 1]);
subplot(2,4,4); bar(T2S_Vector, Max_T2S, 'FaceColor', 'k'); xlabel('T_{2S}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(T2_S, 'r--'); xlim([Lower(7) Upper(7)]); %ylim([0.4 1]);
subplot(2,4,5); bar(kFS_Vector, Max_kFS, 'FaceColor', 'k'); xlabel('k_{FS}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(k_FS, 'r--'); xlim([Lower(5) Upper(5)]); %ylim([0.4 1]);
subplot(2,4,6); bar(kSF_Vector, Max_kSF, 'FaceColor', 'k'); xlabel('k_{SF}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(k_SF, 'r--'); xlim([Lower(6) Upper(6)]); %ylim([0.4 1]);
subplot(2,4,7); bar(M0F_Vector, Max_M0F, 'FaceColor', 'k'); xlabel('M_{0F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(M0_F, 'r--'); xlim([Lower(3) Upper(3)]); %ylim([0.4 1]);
subplot(2,4,8); bar(M0S_Vector, Max_M0S, 'FaceColor', 'k'); xlabel('M_{0S}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); hold on; vline(M0_S, 'r--'); xlim([Lower(4) Upper(4)]); %ylim([0.4 1]);
