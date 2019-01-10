%% Performs MC simulations for histograms.

close all; clear all;

%T1S_1 = 1; T1F_1 = 0.35; T2S_1 = 0.080; T2F_1 = 0.015; k_FS = 0:1:20; M0F_1 = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
%T1S_2 = 1.15; T1F_2 = 0.4; T2S_2 = 0.11; T2F_2 = 0.02; M0F_2 = 0.175; 
%T1S_3 = 1.3; T1F_3 = 0.45; T2S_3 = 0.14; T2F_3 = 0.025; M0F_3 = 0.1; 

T1_S = 1.4; T1_F = 0.45; T2_S = 0.090; T2_F = 0.015; k_FS = 0:1:20; M0_F = 0.15; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

TR_SPGR = 6.5e-3; TR_SSFP = 6.5e-3; FA_SPGR = deg2rad([2 4 6 8 10 12 14 16 18 20]); FA_SSFP180 = deg2rad([2 6 14 22 30 38 46 54 62 70]); FA_SSFP0 = deg2rad([2 6 14 22 30 38 46 54 62 70]);

%% Perform stochastic region contraction.

Realisations = 1; Trials = 40000; Iterations = 30; N = 50; Runs = 1; Params = 5;

Solution_B = zeros(length(k_FS),Realisations,Params);
Solution_D = zeros(length(k_FS),Realisations,Params);
Solution_W = zeros(length(k_FS),Realisations,Params);
Solution_Z = zeros(length(k_FS),Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

parfor kk = 1:length(k_FS)
    
    disp(['k_FS = ', num2str(k_FS(kk))])

    % Ground-truth signals for mcDESPOT.
    SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS(kk));
    
    % Concatenate SSFP signals.
    SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

    for tt = 1:Realisations

        SPGR_Data_Norm = SPGR_Data./mean(SPGR_Data);
        SSFP_Data_Norm = SSFP_Data./mean(SSFP_Data);
        Data_Noiseless = [SPGR_Data_Norm ; SSFP_Data_Norm];

        [Solution_B(kk,tt,:), ~, ~] = SRC_Sim_NoExB(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
        [Solution_D(kk,tt,:), ~, ~] = SRC_Sim_NoExD(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
        [Solution_W(kk,tt,:), ~, ~] = SRC_Sim_NoExW(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
        [Solution_Z(kk,tt,:), ~, ~] = SRC_Sim_NoExZ(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_Noiseless);
                
    end

end

%% Noiseless line plots.

cm = colormap(cool(4));

subplot(2,3,1); plot(k_FS, Solution_B(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 14); ylabel('T_{1S} (s)', 'FontSize', 14); hline(T1_S, 'r--'); grid on; grid minor
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);
subplot(2,3,2); plot(k_FS, Solution_B(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 14); ylabel('T_{1F} (s)', 'FontSize', 14); hline(T1_F, 'r--'); grid on; grid minor
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);
subplot(2,3,3); plot(k_FS, Solution_B(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 14); ylabel('T_{2S} (s)', 'FontSize', 14); hline(T2_S, 'r--'); grid on; grid minor
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);
subplot(2,3,4); plot(k_FS, Solution_B(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 14); ylabel('T_{2F} (s)', 'FontSize', 14); hline(T2_F, 'r--'); grid on; grid minor
get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14);
subplot(2,3,5); plot(k_FS, Solution_B(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 14); ylabel('MWF', 'FontSize', 14); hline(M0_F, 'r--'); grid on; grid minor

subplot(2,3,1); plot(k_FS, Solution_D(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_D(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,3); plot(k_FS, Solution_D(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,4); plot(k_FS, Solution_D(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_D(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  

subplot(2,3,1); plot(k_FS, Solution_W(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_W(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,3); plot(k_FS, Solution_W(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,4); plot(k_FS, Solution_W(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_W(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  

subplot(2,3,1); plot(k_FS, Solution_Z(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,2); plot(k_FS, Solution_Z(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,3); plot(k_FS, Solution_Z(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,4); plot(k_FS, Solution_Z(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,5); plot(k_FS, Solution_Z(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);

% subplot(3,5,6); plot(k_FS, Solution_B_T2(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1S_2, 'r--'); grid on; grid minor
% subplot(3,5,7); plot(k_FS, Solution_B_T2(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1F_2, 'r--'); grid on; grid minor
% subplot(3,5,8); plot(k_FS, Solution_B_T2(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2S_2, 'r--'); grid on; grid minor
% subplot(3,5,9); plot(k_FS, Solution_B_T2(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2F_2, 'r--'); grid on; grid minor
% subplot(3,5,10); plot(k_FS, Solution_B_T2(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0F_2, 'r--'); grid on; grid minor
% 
% subplot(3,5,6); plot(k_FS, Solution_D_T2(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,7); plot(k_FS, Solution_D_T2(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,8); plot(k_FS, Solution_D_T2(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,9); plot(k_FS, Solution_D_T2(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
% subplot(3,5,10); plot(k_FS, Solution_D_T2(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
% 
% subplot(3,5,6); plot(k_FS, Solution_W_T2(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
% subplot(3,5,7); plot(k_FS, Solution_W_T2(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,8); plot(k_FS, Solution_W_T2(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,9); plot(k_FS, Solution_W_T2(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,10); plot(k_FS, Solution_W_T2(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% 
% subplot(3,5,6); plot(k_FS, Solution_Z_T2(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,7); plot(k_FS, Solution_Z_T2(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,8); plot(k_FS, Solution_Z_T2(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,9); plot(k_FS, Solution_Z_T2(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,10); plot(k_FS, Solution_Z_T2(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% 
% subplot(3,5,11); plot(k_FS, Solution_B_T3(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1S_3, 'r--'); grid on; grid minor
% subplot(3,5,12); plot(k_FS, Solution_B_T3(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1F_3, 'r--'); grid on; grid minor
% subplot(3,5,13); plot(k_FS, Solution_B_T3(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2S_3, 'r--'); grid on; grid minor
% subplot(3,5,14); plot(k_FS, Solution_B_T3(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2F_3, 'r--'); grid on; grid minor
% subplot(3,5,15); plot(k_FS, Solution_B_T3(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0F_3, 'r--'); grid on; grid minor
% 
% subplot(3,5,11); plot(k_FS, Solution_D_T3(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,12); plot(k_FS, Solution_D_T3(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,13); plot(k_FS, Solution_D_T3(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
% subplot(3,5,14); plot(k_FS, Solution_D_T3(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
% subplot(3,5,15); plot(k_FS, Solution_D_T3(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
% 
% subplot(3,5,11); plot(k_FS, Solution_W_T3(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
% subplot(3,5,12); plot(k_FS, Solution_W_T3(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,13); plot(k_FS, Solution_W_T3(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,14); plot(k_FS, Solution_W_T3(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% subplot(3,5,15); plot(k_FS, Solution_W_T3(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
% 
% subplot(3,5,11); plot(k_FS, Solution_Z_T3(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,12); plot(k_FS, Solution_Z_T3(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,13); plot(k_FS, Solution_Z_T3(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,14); plot(k_FS, Solution_Z_T3(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
% subplot(3,5,15); plot(k_FS, Solution_Z_T3(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);

ll = legend('B1','B2','B3','B4'); ll.FontSize = 14;
