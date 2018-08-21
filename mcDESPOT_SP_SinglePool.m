%% Performs heat map simulations.

close all; clear all;

T1 = 1; T2 = 0.1; M0 = 1; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

nTrials = 200e6;

P_All = zeros(nTrials,1);
T1_Rand = zeros(nTrials,1); T2_Rand = zeros(nTrials,1); M0_Rand = zeros(nTrials,1);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SP_SteadyState(FA_SPGR, TR_SPGR,'T1',T1,'M0',M0);
SSFP_Data_0 = SSFP_SP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1',T1,'T2',T2,'M0',M0);
SSFP_Data_180 = SSFP_SP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1',T1,'T2',T2,'M0',M0);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

%% Direct sampling of cost-function in 2D.

SNR = 30; Sigma = mean(SPGR_Data)/SNR; % Depends on param set.

Upper = [1.5 1.5 0.15]; Lower = [0.5 0.5 0.05];

% Pre-normalisation.
Data_O = [SPGR_Data ; SSFP_Data];

for nn = 1:nTrials
    
    %disp(['nTrial: ', num2str(nn)])
    
    % Draw random parameter values from a uniform distribution.
    T1_Rand(nn) = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
    M0_Rand(nn) = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
    T2_Rand(nn) = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
    
    P_All(nn) = (logpdf_SinglePool([T1_Rand(nn) M0_Rand(nn) T2_Rand(nn)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
    
end

Exp_P_All = exp(P_All);

%% Signal curve analysis.

PlottingNo = 1000;

[Values, Indexes] = sort(Exp_P_All,'ascend');
Idx_Picked = Indexes(nTrials-(PlottingNo-1):nTrials,1);
T1_Picked = T1_Rand(Idx_Picked,1); T2_Picked = T2_Rand(Idx_Picked,1); M0_Picked = M0_Rand(Idx_Picked,1);

Values_Top = Values(nTrials-(PlottingNo-1):nTrials,1); TSNE_Matrix = [T1_Picked, T2_Picked, M0_Picked];
