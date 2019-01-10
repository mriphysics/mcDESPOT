%%% Performs heat map simulations.

close all; clear all;

acquisition = 'Bouhrara';
exchange = 'on';

switch exchange
    case 'on'
        T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
    case 'off'
        T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 0; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
end

switch acquisition
    case 'Deoni'
        TR_SPGR = 5.6e-3; TR_SSFP = 4.4e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]); FA_SSFP0 = deg2rad([12 16 19 23 27 34 50 70]); FA_SSFP180 = deg2rad([12 16 19 23 27 34 50 70]);
    case 'Bouhrara'
        TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);
end

%% Perform stochastic region contraction.

Realisations = 1; Trials = 40000; Iterations = 30; N = 50; Runs = 1; 

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

SNR = 30; Sigma = mean(SPGR_Data)/SNR;

SPGR_Data_Norm = SPGR_Data./mean(SPGR_Data);
SSFP_Data_Norm = SSFP_Data./mean(SSFP_Data);
Data_Noiseless = [SPGR_Data_Norm ; SSFP_Data_Norm];

switch exchange
    case 'on'
        [Solution, Bounds, Subsets] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
    case 'off'
        [Solution, Bounds, Subsets] = SRC_Sim_NoEx(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
end

%% Direct sampling of cost-function in 2D.

% Pre-normalised data.
Data_O = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Steps = 499;

switch exchange
    
    case 'on'

        Upper = [2 0.7 0.3 16 0.15 0.03]; Lower = [0.8 0.2 0 0 0.03 0];
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps); kFS_Vector = linspace(Lower(4),Upper(4),Steps);
        T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);
        
        P_T1_SRC = zeros(Steps,Steps); P_T2_SRC = zeros(Steps,Steps);
        P_T1F_SRC = zeros(Steps,Steps); P_T1S_SRC = zeros(Steps,Steps);
        P_T2F_SRC = zeros(Steps,Steps); P_T2S_SRC = zeros(Steps,Steps);
        
        for ii = 1:Steps
            for jj = 1:Steps
                
                T1F_Rand = Solution(2); T1S_Rand = Solution(1); T2F_Rand = Solution(4); T2S_Rand = Solution(3); M0F_Rand = Solution(5); kFS_Rand = Solution(6);
                
                P_T1_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Rand T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T1F_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Rand T1F_Vector(ii) M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T1S_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Rand M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2F_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Rand T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2S_SRC(ii,jj) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Vector(jj) T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                
            end
            
        end
        
    case 'off'
        
        Upper = [2 0.7 0.3 0.15 0.03]; Lower = [0.8 0.2 0 0.03 0];
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps);
        T2S_Vector = linspace(Lower(4),Upper(4),Steps); T2F_Vector = linspace(Lower(5),Upper(5),Steps);
        
        P_T1_SRC = zeros(Steps,Steps); P_T2_SRC = zeros(Steps,Steps);
        P_T1F_SRC = zeros(Steps,Steps); P_T1S_SRC = zeros(Steps,Steps);
        P_T2F_SRC = zeros(Steps,Steps); P_T2S_SRC = zeros(Steps,Steps);
        
        for ii = 1:Steps
            for jj = 1:Steps
                
                T1F_Rand = Solution(2); T1S_Rand = Solution(1); T2F_Rand = Solution(4); T2S_Rand = Solution(3); M0F_Rand = Solution(5);
                
                P_T1_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Rand T1F_Rand M0F_Rand T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T1F_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Rand T1F_Vector(ii) M0F_Vector(jj) T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T1S_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Vector(ii) T1F_Rand M0F_Vector(jj) T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2F_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Rand T1F_Rand M0F_Vector(ii) T2S_Rand T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                P_T2S_SRC(ii,jj) = (logpdf_mcDESPOT_NoEx([T1S_Rand T1F_Rand M0F_Vector(ii) T2S_Vector(jj) T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
                
            end
            
        end
        
end
