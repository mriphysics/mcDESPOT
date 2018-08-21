%% Performs heat map simulations.

close all; clear all;

T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

%T1_S = 1; T1_F = 0.35; T2_S = 0.080; T2_F = 0.015; k_FS = 10; M0_F = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta; %10
%T1_S = 1.15; T1_F = 0.4; T2_S = 0.11; T2_F = 0.02; k_FS = 7.5; M0_F = 0.175; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta; %7.5
%T1_S = 1.3; T1_F = 0.45; T2_S = 0.14; T2_F = 0.025; k_FS = 5; M0_F = 0.1; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta; %5

%T1_S = 1; T1_F = 0.35; T2_S = 0.080; T2_F = 0.015; k_FS = 5; M0_F = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
%T1_S = 1; T1_F = 0.35; T2_S = 0.080; T2_F = 0.015; k_FS = 7.5; M0_F = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
%T1_S = 1; T1_F = 0.35; T2_S = 0.080; T2_F = 0.015; k_FS = 10; M0_F = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

%TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]);
%TR_SPGR = 5.6e-3; TR_SSFP = 4.4e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]); FA_SSFP0 = deg2rad([12 16 19 23 27 34 50 70]); FA_SSFP180 = deg2rad([12 16 19 23 27 34 50 70]);
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

nTrials = 200e6;

P_All = zeros(nTrials,1); 
T1S_Rand = zeros(nTrials,1); T1F_Rand = zeros(nTrials,1);
T2S_Rand = zeros(nTrials,1); T2F_Rand = zeros(nTrials,1);
M0F_Rand = zeros(nTrials,1); kFS_Rand = zeros(nTrials,1); Delta_Rand = zeros(nTrials,1);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

%% Direct sampling of cost-function in 2D.

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

SNR = 30; Sigma = mean(SPGR_Data)/SNR;

Upper = [5 0.8 0.5 40 0.2 0.045 (2*pi)]; Lower = [0.7 0.2 0 0.5 0.04 0.001 0]; %Upper = [2 0.7 0.5 20 0.16 0.04]; Lower = [0.8 0.2 0 0.5 0.06 0.002];

% Pre-normalisation.
Data_O = [SPGR_Data ; SSFP_Data];

parfor nn = 1:nTrials
    
    disp(['nTrial: ', num2str(nn)])
    
    % Draw random parameter values from a uniform distribution.
    T1S_Rand(nn) = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
    T1F_Rand(nn) = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
    M0F_Rand(nn) = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
    kFS_Rand(nn) = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
    T2S_Rand(nn) = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
    T2F_Rand(nn) = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);
    Delta_Rand(nn) = (Upper(7) - Lower(7)) .* rand(1,1) + Lower(7);
    
    P_All(nn) = (logpdf_mcDESPOT([T1S_Rand(nn) T1F_Rand(nn) M0F_Rand(nn) kFS_Rand(nn) T2S_Rand(nn) T2F_Rand(nn) Delta_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
    
end

Exp_P_All = exp(P_All);

%% Signal curve analysis.

PlottingNo = 1000;

[Values, Indexes] = sort(Exp_P_All,'ascend');
Idx_Picked = Indexes(nTrials-(PlottingNo-1):nTrials,1);
Values_Picked = Values(nTrials-(PlottingNo-1):nTrials,1);
T1F_Picked = T1F_Rand(Idx_Picked,1); T1S_Picked = T1S_Rand(Idx_Picked,1);
T2F_Picked = T2F_Rand(Idx_Picked,1); T2S_Picked = T2S_Rand(Idx_Picked,1);
M0F_Picked = M0F_Rand(Idx_Picked,1); kFS_Picked = kFS_Rand(Idx_Picked,1); Delta_Picked = Delta_Rand(Idx_Picked,1);
