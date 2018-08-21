%% Performs heat map simulations.

close all; clear all;

% Tissue parameters.
T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta; PC3 = pi/2 + Delta; PC4 = (3*pi)/2 + Delta;

% Acquisition scheme.
TR_SPGR = 5.6e-3; TR_SSFP = 4.4e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]); FA_SSFP0 = deg2rad([12 16 19 23 27 34 50 70]); FA_SSFP180 = deg2rad([12 16 19 23 27 34 50 70]);
%TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]);
%TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

%% Perform stochastic region contraction.

%Realisations = 1; Trials = 40000; Iterations = 30; N = 50; Runs = 1; Params = 6; Solution = zeros(Realisations,Params);
%delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

SNR = 30; Sigma = mean(SPGR_Data)/SNR;

% for tt = 1:Realisations
% 
% %     % Add different noise to each element.
% %     SPGR_Data_Noisy = zeros(length(SPGR_Data),1);
% %     for mm = 1:length(SPGR_Data)
% %         SPGR_Data_Noisy(mm) = SPGR_Data(mm) + (normrnd(0,Sigma));
% %     end
% %     SSFP_Data_Noisy = zeros(length(SSFP_Data),1);
% %     for nn = 1:length(SSFP_Data)
% %         SSFP_Data_Noisy(nn) = SSFP_Data(nn) + (normrnd(0,Sigma));
% %     end
% %     % Normalise signals and concatenate.
% %     SPGR_Data_Norm = SPGR_Data_Noisy./mean(SPGR_Data_Noisy);
% %     SSFP_Data_Norm = SSFP_Data_Noisy./mean(SSFP_Data_Noisy);
% %     Data_NN = [SPGR_Data_Norm ; SSFP_Data_Norm];
%     
%     SPGR_Data_Norm = SPGR_Data./mean(SPGR_Data);
%     SSFP_Data_Norm = SSFP_Data./mean(SSFP_Data);
%     Data_Noiseless = [SPGR_Data_Norm ; SSFP_Data_Norm];
%     
%     % Input post-normalised data.
%     [Solution(tt,:), Bounds, Subsets] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, FA_SSFP90, FA_SSFP270, TR_SPGR, TR_SSFP, Data_Noiseless);
%     %[Solution(tt,:), Bounds, Subsets] = SRC_Sim_NoEx(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
% 
% end

%% Direct sampling of cost-function in 2D. Modify for signal fit plots.

Upper = [2 0.7 0.5 20 0.16 0.04]; Lower = [0.8 0.2 0 0.5 0.06 0.002];

% Pre-normalised data.
Data_O = [SPGR_Data ; SSFP_Data_0 ; SSFP_Data_180];

Steps = 30; nTrials = 50000;

T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); kFS_Vector = linspace(Lower(4),Upper(4),Steps);
T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);
%Delta_Vector = linspace(Lower(7),Upper(7),Steps);

P_T1 = zeros(Steps,Steps,nTrials); P_T2 = zeros(Steps,Steps,nTrials);
P_T1F = zeros(Steps,Steps,nTrials); P_T1S = zeros(Steps,Steps,nTrials);
P_T2F = zeros(Steps,Steps,nTrials); P_T2S = zeros(Steps,Steps,nTrials);
%P_kFS = zeros(Steps,Steps,nTrials);

for ii = 1:Steps
    disp(['Outer Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    for jj = 1:Steps
        %disp(['     Inner Step Number: ', num2str(jj), '/', num2str(Steps), '.'])
        tic
        for nn = 1:nTrials
            
            % Choose GT as first combination.
            %if nn == 1
            %    T1F_Rand = T1_F; T1S_Rand = T1_S; T2F_Rand = T2_F; T2S_Rand = T2_S; kFS_Rand = k_FS; M0F_Rand = M0_F; Delta_Rand = Delta;
            %    %T1F_Rand = Solution(2); T1S_Rand = Solution(1); T2F_Rand = Solution(4); T2S_Rand = Solution(3); M0F_Rand = Solution(5); %kFS_Rand = Solution(6); %Delta_Rand = Solution(7);
            %else
                
                % Draw random parameter values from a uniform distribution.
                T1S_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
                T1F_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
                M0F_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
                kFS_Rand = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
                T2S_Rand = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
                T2F_Rand = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);
                %Delta_Rand = (Upper(7) - Lower(7)) .* rand(1,1) + Lower(7);
                
            %end
            
            P_T1(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma)); 
            P_T2(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Rand kFS_Rand T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1F(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Vector(ii) M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T1S(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Vector(ii) T1F_Rand M0F_Vector(jj) kFS_Rand T2S_Rand T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2F(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Rand T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            P_T2S(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Rand T2S_Vector(jj) T2F_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            %P_kFS(ii,jj,nn) = (logpdf_mcDESPOT([T1S_Rand T1F_Rand M0F_Vector(ii) kFS_Vector(jj) T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_O, Sigma));
            
        end
    end
    
end
