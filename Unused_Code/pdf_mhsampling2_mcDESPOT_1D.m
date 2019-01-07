%% Generates heat maps to depict cost-function space. MaxAngle is in degrees and phase-cycling is in radians.

close all; clear all

% Tissue and sequence parameters.
TR_SPGR = 5e-3; TR_SSFP = 5e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([10 13 17 20 23 30 43 60]);
T1_S = 0.965; T1_F = 0.465; T2_S = 0.09; T2_F = 0.012; M0_F = 0.2; k_FS = 8; M0_S = 0.8; k_SF = (M0_F*k_FS)/M0_S; 
PC1 = 0; PC2 = pi; Sigma = 1/300; %SNR = 30.

%% Ground-truth signals. 

% mcDESPOT GT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SF',k_SF);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SF',k_SF);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SF',k_SF);
Data = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

% Add noise to data.
% Data_Noisy = zeros(length(Data),1); Sigma_Data = 1/300;
% for jj = 1:length(Data)
%     Data_Noisy(jj) = Data(jj) + (normrnd(0,Sigma_Data));
% end

Upper = [1.25 0.5 0.30 20 0.15 0.03]; Lower = [0.20 0.1 0.05 1 0.04 0.01];

%% Direct sampling of cost-function in 2D.

Steps = 100; nTrials = 20000;
T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); kFS_Vector = linspace(Lower(4),Upper(4),Steps);
T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);

P_T1F = zeros(length(T1F_Vector),1); P_T1S = zeros(length(T1S_Vector),1);
P_kFS = zeros(length(kFS_Vector),1); P_M0F = zeros(length(M0F_Vector),1);
P_T2F = zeros(length(T2F_Vector),1); P_T2S = zeros(length(T2S_Vector),1);

for ii = 1:length(T1S_Vector)
    tic
    disp(['Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    %%% CHANGE INPUT DATA TO DATA_NOISY? %%%
    for nn = 1:nTrials
        
        % Draw random parameter values from a uniform distribution.
        T1S_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
        T1F_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
        M0F_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
        kFS_Rand = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
        T2S_Rand = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
        T2F_Rand = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);

        P_T1F(ii,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Vector(ii) M0F_Rand kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1S(ii,nn) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Rand M0F_Rand kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Vector(ii) T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_M0F(ii,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Rand T2S_Rand T2F_Vector(ii)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Rand T2S_Vector(ii) T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
   
    end
    toc
end

Max_T1S = max(exp(P_T1S),[],2); Max_T1F = max(exp(P_T1F),[],2);
Max_T2S = max(exp(P_T2S),[],2); Max_T2F = max(exp(P_T2F),[],2);
Max_kFS = max(exp(P_kFS),[],2); Max_M0F = max(exp(P_M0F),[],2); 

subplot(2,3,1); bar(T1F_Vector, Max_T1F, 'FaceColor', 'k'); xlabel('T_{1F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(T1_F, 'r--'); xlim([Lower(2) Upper(2)]); ylim([0.4 1]);
subplot(2,3,2); bar(T1S_Vector, Max_T1S, 'FaceColor', 'k'); xlabel('T_{1S}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(T1_S, 'r--'); xlim([Lower(1) Upper(1)]); ylim([0.4 1]);
subplot(2,3,3); bar(T2F_Vector, Max_T2F, 'FaceColor', 'k'); xlabel('T_{2F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(T2_F, 'r--'); xlim([Lower(6) Upper(6)]); ylim([0.4 1]);
subplot(2,3,4); bar(T2S_Vector, Max_T2S, 'FaceColor', 'k'); xlabel('T_{2S}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(T2_S, 'r--'); xlim([Lower(5) Upper(5)]); ylim([0.4 1]);
subplot(2,3,5); bar(kFS_Vector, Max_kFS, 'FaceColor', 'k'); xlabel('k_{FS}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(k_FS, 'r--'); xlim([Lower(4) Upper(4)]); ylim([0.4 1]);
subplot(2,3,6); bar(M0F_Vector, Max_M0F, 'FaceColor', 'k'); xlabel('M_{0F}','FontSize',12); ylabel('Max. Likelihood', 'FontSize', 12); vline(M0_F, 'r--'); xlim([Lower(3) Upper(3)]); ylim([0.4 1]);
