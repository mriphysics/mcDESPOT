%% Performs MC simulations for histograms.

%close all; clear all;

T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; k_FS = 0:1:20; M0F_1 = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; M0F_2 = 0.175; 
T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; M0F_3 = 0.1; 

TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

%% Perform stochastic region contraction.

Realisations = 1; Trials = 40000; Iterations = 30; N = 50; Runs = 1; Params = 5;

Solution = zeros(length(k_FS),Realisations,Params);

%delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

for kk = 1:length(k_FS)
    
    disp(['k_FS = ', num2str(k_FS(kk))])

    % Ground-truth signals for mcDESPOT.
    SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    
    % Concatenate SSFP signals.
    SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];
    
    % Define SNR and define wrt mean SPGR signal.
    %SNR = 30; Sigma = mean(SPGR_Data)/SNR;
    
    for tt = 1:Realisations
        
        %disp(['Realisation: ', num2str(tt)])
        
        % Add different noise to each element.
        %SPGR_Data_Noisy = zeros(length(SPGR_Data),1);
        %for mm = 1:length(SPGR_Data)
        %   SPGR_Data_Noisy(mm) = SPGR_Data(mm) + (normrnd(0,Sigma));
        %end
        %SSFP_Data_Noisy = zeros(length(SSFP_Data),1);
        %for nn = 1:length(SSFP_Data)
        %   SSFP_Data_Noisy(nn) = SSFP_Data(nn) + (normrnd(0,Sigma));
        %end
        
        % Normalise signals and concatenate.
        %SPGR_Data_Norm = SPGR_Data_Noisy./mean(SPGR_Data_Noisy);
        %SSFP_Data_Norm = SSFP_Data_Noisy./mean(SSFP_Data_Noisy);
        %Data_NN = [SPGR_Data_Norm ; SSFP_Data_Norm];
        
        SPGR_Data_Norm = SPGR_Data./mean(SPGR_Data);
        SSFP_Data_Norm = SSFP_Data./mean(SSFP_Data);
        Data_Noiseless = [SPGR_Data_Norm ; SSFP_Data_Norm];
        
        % Post-normalisation.
        %[Solution(tt,:), ~, ~] = SRC_Sim(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_NN);
        [Solution(kk,tt,:), ~, ~] = SRC_Sim_NoEx_3(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
        
    end

end

%% Noiseless line plots.

cm = colormap(cool(4));

subplot(3,5,1); plot(k_FS, Solution_B_T1(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1S_1, 'r--'); grid on; grid minor
subplot(3,5,2); plot(k_FS, Solution_B_T1(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1F_1, 'r--'); grid on; grid minor
subplot(3,5,3); plot(k_FS, Solution_B_T1(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2S_1, 'r--'); grid on; grid minor
subplot(3,5,4); plot(k_FS, Solution_B_T1(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2F_1, 'r--'); grid on; grid minor
subplot(3,5,5); plot(k_FS, Solution_B_T1(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0F_1, 'r--'); grid on; grid minor

subplot(3,5,1); plot(k_FS, Solution_D_T1(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,2); plot(k_FS, Solution_D_T1(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,3); plot(k_FS, Solution_D_T1(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,4); plot(k_FS, Solution_D_T1(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(3,5,5); plot(k_FS, Solution_D_T1(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  

subplot(3,5,1); plot(k_FS, Solution_W_T1(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(3,5,2); plot(k_FS, Solution_W_T1(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,3); plot(k_FS, Solution_W_T1(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,4); plot(k_FS, Solution_W_T1(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,5); plot(k_FS, Solution_W_T1(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  

subplot(3,5,1); plot(k_FS, Solution_Z_T1(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,2); plot(k_FS, Solution_Z_T1(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,3); plot(k_FS, Solution_Z_T1(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,4); plot(k_FS, Solution_Z_T1(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,5); plot(k_FS, Solution_Z_T1(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);

subplot(3,5,6); plot(k_FS, Solution_B_T2(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1S_2, 'r--'); grid on; grid minor
subplot(3,5,7); plot(k_FS, Solution_B_T2(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1F_2, 'r--'); grid on; grid minor
subplot(3,5,8); plot(k_FS, Solution_B_T2(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2S_2, 'r--'); grid on; grid minor
subplot(3,5,9); plot(k_FS, Solution_B_T2(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2F_2, 'r--'); grid on; grid minor
subplot(3,5,10); plot(k_FS, Solution_B_T2(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0F_2, 'r--'); grid on; grid minor

subplot(3,5,6); plot(k_FS, Solution_D_T2(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,7); plot(k_FS, Solution_D_T2(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,8); plot(k_FS, Solution_D_T2(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,9); plot(k_FS, Solution_D_T2(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(3,5,10); plot(k_FS, Solution_D_T2(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  

subplot(3,5,6); plot(k_FS, Solution_W_T2(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(3,5,7); plot(k_FS, Solution_W_T2(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,8); plot(k_FS, Solution_W_T2(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,9); plot(k_FS, Solution_W_T2(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,10); plot(k_FS, Solution_W_T2(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  

subplot(3,5,6); plot(k_FS, Solution_Z_T2(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,7); plot(k_FS, Solution_Z_T2(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,8); plot(k_FS, Solution_Z_T2(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,9); plot(k_FS, Solution_Z_T2(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,10); plot(k_FS, Solution_Z_T2(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);

subplot(3,5,11); plot(k_FS, Solution_B_T3(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1S_3, 'r--'); grid on; grid minor
subplot(3,5,12); plot(k_FS, Solution_B_T3(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1F_3, 'r--'); grid on; grid minor
subplot(3,5,13); plot(k_FS, Solution_B_T3(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2S_3, 'r--'); grid on; grid minor
subplot(3,5,14); plot(k_FS, Solution_B_T3(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2F_3, 'r--'); grid on; grid minor
subplot(3,5,15); plot(k_FS, Solution_B_T3(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0F_3, 'r--'); grid on; grid minor

subplot(3,5,11); plot(k_FS, Solution_D_T3(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,12); plot(k_FS, Solution_D_T3(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,13); plot(k_FS, Solution_D_T3(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(3,5,14); plot(k_FS, Solution_D_T3(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(3,5,15); plot(k_FS, Solution_D_T3(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  

subplot(3,5,11); plot(k_FS, Solution_W_T3(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(3,5,12); plot(k_FS, Solution_W_T3(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,13); plot(k_FS, Solution_W_T3(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,14); plot(k_FS, Solution_W_T3(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(3,5,15); plot(k_FS, Solution_W_T3(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  

subplot(3,5,11); plot(k_FS, Solution_Z_T3(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,12); plot(k_FS, Solution_Z_T3(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,13); plot(k_FS, Solution_Z_T3(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,14); plot(k_FS, Solution_Z_T3(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(3,5,15); plot(k_FS, Solution_Z_T3(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);

ll = legend('Set 1','Set 2','Set 3','Set 4'); ll.FontSize = 12;
