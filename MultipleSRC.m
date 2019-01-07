
function Solution = MultipleSRC(T1_F, T1_S, T2_F, T2_S, M0_F, k_FS)

TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

%% Perform stochastic region contraction.

Realisations = 1000; Trials = 40000; Iterations = 30; N = 50; Runs = 1; Params = 6;

Solution = zeros(Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

% Define SNR and define wrt mean SPGR signal.
SNR = 30; Sigma = mean(SPGR_Data)/SNR;

parfor tt = 1:Realisations
        
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
    Data_NN = [SPGR_Data_Norm ; SSFP_Data_Norm];

    % Post-normalisation.
    [Solution(tt,:), ~, ~] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_NN);
    
end

end