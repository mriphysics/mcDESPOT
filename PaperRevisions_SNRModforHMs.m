%%% SNR = 30 --> SNR = 100 conversion script for search-space results. %%%

%% Exchange, Figure 2.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 8; nMeas = 30; 
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

E_Values_Top_SNR100 = exp((log(E_Values_Top)*SNR_Factor));
E_Values_Top_NoiseInd = sqrt((-2*log(E_Values_Top_SNR100))/nMeas);

%% No exchange, Figure 2.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 0; nMeas = 30; 
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

NE_Values_Top_SNR100 = exp((log(NE_Values_Top)*SNR_Factor));
NE_Values_Top_NoiseInd = sqrt((-2*log(NE_Values_Top_SNR100))/nMeas);

%% Single-pool, Figure 2.

T1 = 1; T2 = 0.1; M0 = 1; nMeas = 30; 
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

SP_Values_Top_SNR100 = exp((log(SP_Values_Top)*SNR_Factor));
SP_Values_Top_NoiseInd = sqrt((-2*log(SP_Values_Top_SNR100))/nMeas);

%% Exchange, Figure 3.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 8; nMeas = 30; 
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

MaxCF_T1_GT_SNR100 = exp((log(MaxCF_T1_GT)*SNR_Factor));
MaxCF_T1_GT_NoiseInd = sqrt((-2*log(MaxCF_T1_GT_SNR100))/nMeas);
MaxCF_T1_MP_SNR100 = exp((log(MaxCF_T1_MP)*SNR_Factor));
MaxCF_T1_MP_NoiseInd = sqrt((-2*log(MaxCF_T1_MP_SNR100))/nMeas);
MaxCF_T1_SRC_SNR100 = exp((log(MaxCF_T1_SRC)*SNR_Factor));
MaxCF_T1_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1_SRC_SNR100))/nMeas);

MaxCF_T1F_GT_SNR100 = exp((log(MaxCF_T1F_GT)*SNR_Factor));
MaxCF_T1F_GT_NoiseInd = sqrt((-2*log(MaxCF_T1F_GT_SNR100))/nMeas);
MaxCF_T1F_MP_SNR100 = exp((log(MaxCF_T1F_MP)*SNR_Factor));
MaxCF_T1F_MP_NoiseInd = sqrt((-2*log(MaxCF_T1F_MP_SNR100))/nMeas);
MaxCF_T1F_SRC_SNR100 = exp((log(MaxCF_T1F_SRC)*SNR_Factor));
MaxCF_T1F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1F_SRC_SNR100))/nMeas);

MaxCF_T1S_GT_SNR100 = exp((log(MaxCF_T1S_GT)*SNR_Factor));
MaxCF_T1S_GT_NoiseInd = sqrt((-2*log(MaxCF_T1S_GT_SNR100))/nMeas);
MaxCF_T1S_MP_SNR100 = exp((log(MaxCF_T1S_MP)*SNR_Factor));
MaxCF_T1S_MP_NoiseInd = sqrt((-2*log(MaxCF_T1S_MP_SNR100))/nMeas);
MaxCF_T1S_SRC_SNR100 = exp((log(MaxCF_T1S_SRC)*SNR_Factor));
MaxCF_T1S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1S_SRC_SNR100))/nMeas);

MaxCF_T2_GT_SNR100 = exp((log(MaxCF_T2_GT)*SNR_Factor));
MaxCF_T2_GT_NoiseInd = sqrt((-2*log(MaxCF_T2_GT_SNR100))/nMeas);
MaxCF_T2_MP_SNR100 = exp((log(MaxCF_T2_MP)*SNR_Factor));
MaxCF_T2_MP_NoiseInd = sqrt((-2*log(MaxCF_T2_MP_SNR100))/nMeas);
MaxCF_T2_SRC_SNR100 = exp((log(MaxCF_T2_SRC)*SNR_Factor));
MaxCF_T2_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2_SRC_SNR100))/nMeas);

MaxCF_T2F_GT_SNR100 = exp((log(MaxCF_T2F_GT)*SNR_Factor));
MaxCF_T2F_GT_NoiseInd = sqrt((-2*log(MaxCF_T2F_GT_SNR100))/nMeas);
MaxCF_T2F_MP_SNR100 = exp((log(MaxCF_T2F_MP)*SNR_Factor));
MaxCF_T2F_MP_NoiseInd = sqrt((-2*log(MaxCF_T2F_MP_SNR100))/nMeas);
MaxCF_T2F_SRC_SNR100 = exp((log(MaxCF_T2F_SRC)*SNR_Factor));
MaxCF_T2F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2F_SRC_SNR100))/nMeas);

MaxCF_T2S_GT_SNR100 = exp((log(MaxCF_T2S_GT)*SNR_Factor));
MaxCF_T2S_GT_NoiseInd = sqrt((-2*log(MaxCF_T2S_GT_SNR100))/nMeas);
MaxCF_T2S_MP_SNR100 = exp((log(MaxCF_T2S_MP)*SNR_Factor));
MaxCF_T2S_MP_NoiseInd = sqrt((-2*log(MaxCF_T2S_MP_SNR100))/nMeas);
MaxCF_T2S_SRC_SNR100 = exp((log(MaxCF_T2S_SRC)*SNR_Factor));
MaxCF_T2S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2S_SRC_SNR100))/nMeas);


%% No exchange, Figure 4.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 0; nMeas = 30;
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

MaxCF_T1_GT_SNR100 = exp((log(MaxCF_T1_GT)*SNR_Factor));
MaxCF_T1_GT_NoiseInd = sqrt((-2*log(MaxCF_T1_GT_SNR100))/nMeas);
MaxCF_T1_MP_SNR100 = exp((log(MaxCF_T1_MP)*SNR_Factor));
MaxCF_T1_MP_NoiseInd = sqrt((-2*log(MaxCF_T1_MP_SNR100))/nMeas);
MaxCF_T1_SRC_SNR100 = exp((log(MaxCF_T1_SRC)*SNR_Factor));
MaxCF_T1_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1_SRC_SNR100))/nMeas);

MaxCF_T1F_GT_SNR100 = exp((log(MaxCF_T1F_GT)*SNR_Factor));
MaxCF_T1F_GT_NoiseInd = sqrt((-2*log(MaxCF_T1F_GT_SNR100))/nMeas);
MaxCF_T1F_MP_SNR100 = exp((log(MaxCF_T1F_MP)*SNR_Factor));
MaxCF_T1F_MP_NoiseInd = sqrt((-2*log(MaxCF_T1F_MP_SNR100))/nMeas);
MaxCF_T1F_SRC_SNR100 = exp((log(MaxCF_T1F_SRC)*SNR_Factor));
MaxCF_T1F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1F_SRC_SNR100))/nMeas);

MaxCF_T1S_GT_SNR100 = exp((log(MaxCF_T1S_GT)*SNR_Factor));
MaxCF_T1S_GT_NoiseInd = sqrt((-2*log(MaxCF_T1S_GT_SNR100))/nMeas);
MaxCF_T1S_MP_SNR100 = exp((log(MaxCF_T1S_MP)*SNR_Factor));
MaxCF_T1S_MP_NoiseInd = sqrt((-2*log(MaxCF_T1S_MP_SNR100))/nMeas);
MaxCF_T1S_SRC_SNR100 = exp((log(MaxCF_T1S_SRC)*SNR_Factor));
MaxCF_T1S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1S_SRC_SNR100))/nMeas);

MaxCF_T2_GT_SNR100 = exp((log(MaxCF_T2_GT)*SNR_Factor));
MaxCF_T2_GT_NoiseInd = sqrt((-2*log(MaxCF_T2_GT_SNR100))/nMeas);
MaxCF_T2_MP_SNR100 = exp((log(MaxCF_T2_MP)*SNR_Factor));
MaxCF_T2_MP_NoiseInd = sqrt((-2*log(MaxCF_T2_MP_SNR100))/nMeas);
MaxCF_T2_SRC_SNR100 = exp((log(MaxCF_T2_SRC)*SNR_Factor));
MaxCF_T2_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2_SRC_SNR100))/nMeas);

MaxCF_T2F_GT_SNR100 = exp((log(MaxCF_T2F_GT)*SNR_Factor));
MaxCF_T2F_GT_NoiseInd = sqrt((-2*log(MaxCF_T2F_GT_SNR100))/nMeas);
MaxCF_T2F_MP_SNR100 = exp((log(MaxCF_T2F_MP)*SNR_Factor));
MaxCF_T2F_MP_NoiseInd = sqrt((-2*log(MaxCF_T2F_MP_SNR100))/nMeas);
MaxCF_T2F_SRC_SNR100 = exp((log(MaxCF_T2F_SRC)*SNR_Factor));
MaxCF_T2F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2F_SRC_SNR100))/nMeas);

MaxCF_T2S_GT_SNR100 = exp((log(MaxCF_T2S_GT)*SNR_Factor));
MaxCF_T2S_GT_NoiseInd = sqrt((-2*log(MaxCF_T2S_GT_SNR100))/nMeas);
MaxCF_T2S_MP_SNR100 = exp((log(MaxCF_T2S_MP)*SNR_Factor));
MaxCF_T2S_MP_NoiseInd = sqrt((-2*log(MaxCF_T2S_MP_SNR100))/nMeas);
MaxCF_T2S_SRC_SNR100 = exp((log(MaxCF_T2S_SRC)*SNR_Factor));
MaxCF_T2S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2S_SRC_SNR100))/nMeas);

%% Exchange, SF 3.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 8; nMeas = 24; 
TR_SPGR = 5.6e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

MaxCF_T1_GT_SNR100 = exp((log(MaxCF_T1_GT)*SNR_Factor));
MaxCF_T1_GT_NoiseInd = sqrt((-2*log(MaxCF_T1_GT_SNR100))/nMeas);
MaxCF_T1_MP_SNR100 = exp((log(MaxCF_T1_MP)*SNR_Factor));
MaxCF_T1_MP_NoiseInd = sqrt((-2*log(MaxCF_T1_MP_SNR100))/nMeas);
MaxCF_T1_SRC_SNR100 = exp((log(MaxCF_T1_SRC)*SNR_Factor));
MaxCF_T1_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1_SRC_SNR100))/nMeas);

MaxCF_T1F_GT_SNR100 = exp((log(MaxCF_T1F_GT)*SNR_Factor));
MaxCF_T1F_GT_NoiseInd = sqrt((-2*log(MaxCF_T1F_GT_SNR100))/nMeas);
MaxCF_T1F_MP_SNR100 = exp((log(MaxCF_T1F_MP)*SNR_Factor));
MaxCF_T1F_MP_NoiseInd = sqrt((-2*log(MaxCF_T1F_MP_SNR100))/nMeas);
MaxCF_T1F_SRC_SNR100 = exp((log(MaxCF_T1F_SRC)*SNR_Factor));
MaxCF_T1F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1F_SRC_SNR100))/nMeas);

MaxCF_T1S_GT_SNR100 = exp((log(MaxCF_T1S_GT)*SNR_Factor));
MaxCF_T1S_GT_NoiseInd = sqrt((-2*log(MaxCF_T1S_GT_SNR100))/nMeas);
MaxCF_T1S_MP_SNR100 = exp((log(MaxCF_T1S_MP)*SNR_Factor));
MaxCF_T1S_MP_NoiseInd = sqrt((-2*log(MaxCF_T1S_MP_SNR100))/nMeas);
MaxCF_T1S_SRC_SNR100 = exp((log(MaxCF_T1S_SRC)*SNR_Factor));
MaxCF_T1S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1S_SRC_SNR100))/nMeas);

MaxCF_T2_GT_SNR100 = exp((log(MaxCF_T2_GT)*SNR_Factor));
MaxCF_T2_GT_NoiseInd = sqrt((-2*log(MaxCF_T2_GT_SNR100))/nMeas);
MaxCF_T2_MP_SNR100 = exp((log(MaxCF_T2_MP)*SNR_Factor));
MaxCF_T2_MP_NoiseInd = sqrt((-2*log(MaxCF_T2_MP_SNR100))/nMeas);
MaxCF_T2_SRC_SNR100 = exp((log(MaxCF_T2_SRC)*SNR_Factor));
MaxCF_T2_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2_SRC_SNR100))/nMeas);

MaxCF_T2F_GT_SNR100 = exp((log(MaxCF_T2F_GT)*SNR_Factor));
MaxCF_T2F_GT_NoiseInd = sqrt((-2*log(MaxCF_T2F_GT_SNR100))/nMeas);
MaxCF_T2F_MP_SNR100 = exp((log(MaxCF_T2F_MP)*SNR_Factor));
MaxCF_T2F_MP_NoiseInd = sqrt((-2*log(MaxCF_T2F_MP_SNR100))/nMeas);
MaxCF_T2F_SRC_SNR100 = exp((log(MaxCF_T2F_SRC)*SNR_Factor));
MaxCF_T2F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2F_SRC_SNR100))/nMeas);

MaxCF_T2S_GT_SNR100 = exp((log(MaxCF_T2S_GT)*SNR_Factor));
MaxCF_T2S_GT_NoiseInd = sqrt((-2*log(MaxCF_T2S_GT_SNR100))/nMeas);
MaxCF_T2S_MP_SNR100 = exp((log(MaxCF_T2S_MP)*SNR_Factor));
MaxCF_T2S_MP_NoiseInd = sqrt((-2*log(MaxCF_T2S_MP_SNR100))/nMeas);
MaxCF_T2S_SRC_SNR100 = exp((log(MaxCF_T2S_SRC)*SNR_Factor));
MaxCF_T2S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2S_SRC_SNR100))/nMeas);


%% No Exchange, SF 3.

T1_F = 0.45; T1_S = 1.4; M0_F = 0.15; k_FS = 0; nMeas = 24;
TR_SPGR = 5.6e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]);
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

MaxCF_T1_GT_SNR100 = exp((log(MaxCF_T1_GT)*SNR_Factor));
MaxCF_T1_GT_NoiseInd = sqrt((-2*log(MaxCF_T1_GT_SNR100))/nMeas);
MaxCF_T1_MP_SNR100 = exp((log(MaxCF_T1_MP)*SNR_Factor));
MaxCF_T1_MP_NoiseInd = sqrt((-2*log(MaxCF_T1_MP_SNR100))/nMeas);
MaxCF_T1_SRC_SNR100 = exp((log(MaxCF_T1_SRC)*SNR_Factor));
MaxCF_T1_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1_SRC_SNR100))/nMeas);

MaxCF_T1F_GT_SNR100 = exp((log(MaxCF_T1F_GT)*SNR_Factor));
MaxCF_T1F_GT_NoiseInd = sqrt((-2*log(MaxCF_T1F_GT_SNR100))/nMeas);
MaxCF_T1F_MP_SNR100 = exp((log(MaxCF_T1F_MP)*SNR_Factor));
MaxCF_T1F_MP_NoiseInd = sqrt((-2*log(MaxCF_T1F_MP_SNR100))/nMeas);
MaxCF_T1F_SRC_SNR100 = exp((log(MaxCF_T1F_SRC)*SNR_Factor));
MaxCF_T1F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1F_SRC_SNR100))/nMeas);

MaxCF_T1S_GT_SNR100 = exp((log(MaxCF_T1S_GT)*SNR_Factor));
MaxCF_T1S_GT_NoiseInd = sqrt((-2*log(MaxCF_T1S_GT_SNR100))/nMeas);
MaxCF_T1S_MP_SNR100 = exp((log(MaxCF_T1S_MP)*SNR_Factor));
MaxCF_T1S_MP_NoiseInd = sqrt((-2*log(MaxCF_T1S_MP_SNR100))/nMeas);
MaxCF_T1S_SRC_SNR100 = exp((log(MaxCF_T1S_SRC)*SNR_Factor));
MaxCF_T1S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T1S_SRC_SNR100))/nMeas);

MaxCF_T2_GT_SNR100 = exp((log(MaxCF_T2_GT)*SNR_Factor));
MaxCF_T2_GT_NoiseInd = sqrt((-2*log(MaxCF_T2_GT_SNR100))/nMeas);
MaxCF_T2_MP_SNR100 = exp((log(MaxCF_T2_MP)*SNR_Factor));
MaxCF_T2_MP_NoiseInd = sqrt((-2*log(MaxCF_T2_MP_SNR100))/nMeas);
MaxCF_T2_SRC_SNR100 = exp((log(MaxCF_T2_SRC)*SNR_Factor));
MaxCF_T2_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2_SRC_SNR100))/nMeas);

MaxCF_T2F_GT_SNR100 = exp((log(MaxCF_T2F_GT)*SNR_Factor));
MaxCF_T2F_GT_NoiseInd = sqrt((-2*log(MaxCF_T2F_GT_SNR100))/nMeas);
MaxCF_T2F_MP_SNR100 = exp((log(MaxCF_T2F_MP)*SNR_Factor));
MaxCF_T2F_MP_NoiseInd = sqrt((-2*log(MaxCF_T2F_MP_SNR100))/nMeas);
MaxCF_T2F_SRC_SNR100 = exp((log(MaxCF_T2F_SRC)*SNR_Factor));
MaxCF_T2F_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2F_SRC_SNR100))/nMeas);

MaxCF_T2S_GT_SNR100 = exp((log(MaxCF_T2S_GT)*SNR_Factor));
MaxCF_T2S_GT_NoiseInd = sqrt((-2*log(MaxCF_T2S_GT_SNR100))/nMeas);
MaxCF_T2S_MP_SNR100 = exp((log(MaxCF_T2S_MP)*SNR_Factor));
MaxCF_T2S_MP_NoiseInd = sqrt((-2*log(MaxCF_T2S_MP_SNR100))/nMeas);
MaxCF_T2S_SRC_SNR100 = exp((log(MaxCF_T2S_SRC)*SNR_Factor));
MaxCF_T2S_SRC_NoiseInd = sqrt((-2*log(MaxCF_T2S_SRC_SNR100))/nMeas);

%% Single-Pool, SF 4.

T1 = 1; T2 = 0.1; M0 = 1; Delta = 0; nMeas = 30;
TR_SPGR = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]);
SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
SNR_1 = 30; Sigma_1 = mean(SPGR_Data)/SNR_1;

SNR_2 = 100; Sigma_2 = mean(SPGR_Data)/SNR_2;

SNR_Factor = (Sigma_1/Sigma_2)^2;

MaxCF_T1_SNR100 = exp((log(MaxCF_T1)*SNR_Factor));
MaxCF_T1_NoiseInd = sqrt((-2*log(MaxCF_T1_SNR100))/nMeas);

MaxCF_T2_SNR100 = exp((log(MaxCF_T2)*SNR_Factor));
MaxCF_T2_NoiseInd = sqrt((-2*log(MaxCF_T2_SNR100))/nMeas);

MaxCF_T1T2_SNR100 = exp((log(MaxCF_T1T2)*SNR_Factor));
MaxCF_T1T2_NoiseInd = sqrt((-2*log(MaxCF_T1T2_SNR100))/nMeas);