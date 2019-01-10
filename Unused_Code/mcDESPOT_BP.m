%% Performs MC simulations for box plots.

close all; clear all;

% Tissue and sequence parameters.
T1_S = 0.965; T1S_Fix = [0.9, 0.965, 1.1, 1.3, 1.4, 1.5];
T1_F = 0.465; T1F_Fix = [0.3, 0.4, 0.465, 0.6, 0.7, 0.8]; 
T2_S = 0.090; T2S_Fix = [0.04, 0.06, 0.09, 0.11, 0.13, 0.15]; 
T2_F = 0.012; T2F_Fix = [0.01, 0.012, 0.015, 0.02, 0.025, 0.03];
k_FS = 8.000; kFS_Fix = [(1/0.6), 5, 8, 10, 20, 40];
M0_F = 0.200; 
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;
% T1_S = 1.2; T1_F = 0.55; T2_S = 0.095; T2_F = 0.02; k_FS = 10; M0_F = 0.175; Delta = pi/8; PC1 = 0 + Delta; PC2 = pi + Delta;
TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]);
%TR_SPGR = 5.6e-3; TR_SSFP = 4.4e-3; FA_SPGR = deg2rad([4 5 6 7 9 11 14 18]); FA_SSFP0 = deg2rad([12 16 19 23 27 34 50 70]); FA_SSFP180 = deg2rad([12 16 19 23 27 34 50 70]); 

% Ground-truth signals for mcDESPOT.
SPGR_Data = SPGR_SteadyState(FA_SPGR, TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_0 = SSFP_SteadyState(FA_SSFP0, TR_SSFP, PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_SteadyState(FA_SSFP180, TR_SSFP, PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);

% Concatenate SSFP signals.
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];

%% Perform stochastic region contraction.

Realisations = 100; Trials = 5000; Iterations = 7; N = 50; Runs = 1; Params = 7;
SNR = 30; Sigma = mean(SPGR_Data)/SNR;

Solution = zeros(1,Realisations,Params);

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 8; parpool(c, c.NumWorkers);

tic
for ii = 1:1

    parfor tt = 1:Realisations
        disp(['Realisation: ', num2str(tt)])
        
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
        [Solution(ii,tt,:), ~, ~] = SRC_Sim_Fix(Trials, Iterations, N, Runs, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data_NN);
        
    end

end
toc

%% MWF box plots.

% BP_Vector = [M0F_All ; M0F_FixedT1F ; M0F_FixedT1S ; M0F_FixedT2F ; M0F_FixedT2S ; M0F_FixedkFS ; M0F_FixedWT1F ; M0F_FixedWT1S ; M0F_FixedWT2F ; M0F_FixedWT2S ; M0F_FixedWkFS ; M0F_FixedT1FT2F ; M0F_FixedT1ST2S ; M0F_FixedT1FT1S ; M0F_FixedT2FT2S];
% Plotting_Vector = [ones(100,1); 2*ones(100,1); 3*ones(100,1); 4*ones(100,1); 5*ones(100,1); 6*ones(100,1); 7*ones(100,1); 8*ones(100,1); 9*ones(100,1); 10*ones(100,1); 11*ones(100,1); 12*ones(100,1); 13*ones(100,1); 14*ones(100,1); 15*ones(100,1)];
% boxplot(BP_Vector,Plotting_Vector)
% Labels = {'All','T1F','T1S','T2F','T2S','kFS','2T1F','2T1S','2T2F','2T2S','2kFS','T1F&T2F','T1S&T2S','T1F&T1S','T2F&T2S'};
% set(gca,'XTickLabel',Labels)
% hold on; hline(0.2, 'm--');

% For Wood's bounds.
% figure(1);
% subplot(2,3,1); boxplot([Solution_All(1,:,5)', Solution_T1F(1,:,5)', Solution_T1F(2,:,5)', Solution_T1F(3,:,5)', Solution_T1F(4,:,5)', Solution_T1F(5,:,5)', Solution_T1F(6,:,5)']);
% a = get(get(gca,'children'),'children'); t = get(a,'tag'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.3','0.4','0.465','0.6','0.7','0.8'}; set(gca,'XTickLabel',Labels); xlabel('T_{1F} (s)'); hline(0.2, 'm--');
% subplot(2,3,2); boxplot([Solution_All(1,:,5)', Solution_T1S(1,:,5)', Solution_T1S(2,:,5)', Solution_T1S(3,:,5)', Solution_T1S(4,:,5)', Solution_T1S(5,:,5)', Solution_T1S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.9','0.965','1.1','1.3','1.4','1.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1S} (s)'); hline(0.2, 'm--');
% subplot(2,3,3); boxplot([Solution_All(1,:,5)', Solution_T2F(1,:,5)', Solution_T2F(2,:,5)', Solution_T2F(3,:,5)', Solution_T2F(4,:,5)', Solution_T2F(5,:,5)', Solution_T2F(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','10','12','15','20','25','30'}; set(gca,'XTickLabel',Labels); xlabel('T_{2F} (ms)'); hline(0.2, 'm--');
% subplot(2,3,4); boxplot([Solution_All(1,:,5)', Solution_T2S(1,:,5)', Solution_T2S(2,:,5)', Solution_T2S(3,:,5)', Solution_T2S(4,:,5)', Solution_T2S(5,:,5)', Solution_T2S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','40','60','90','110','130','150'}; set(gca,'XTickLabel',Labels); xlabel('T_{2S} (ms)'); hline(0.2, 'm--');
% subplot(2,3,5); boxplot([Solution_All(1,:,5)', Solution_kFS(1,:,5)', Solution_kFS(2,:,5)', Solution_kFS(3,:,5)', Solution_kFS(4,:,5)', Solution_kFS(5,:,5)', Solution_kFS(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','(1/0.6)','5','8','10','20','40'}; set(gca,'XTickLabel',Labels); xlabel('k_{FS} (s^{-1})'); hline(0.2, 'm--');
 
% % For Zhang's bounds.
% figure(2);
% subplot(2,3,1); boxplot([Solution(:,5), SolutionT1F(1,:,5)', SolutionT1F(2,:,5)', SolutionT1F(3,:,5)', SolutionT1F(4,:,5)', SolutionT1F(5,:,5)', SolutionT1F(6,:,5)']);
% a = get(get(gca,'children'),'children'); box = a(16); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.2','0.25','0.3','0.4','0.465','0.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1F} (s)'); hline(0.2, 'm--');
% subplot(2,3,2); boxplot([Solution(:,5), SolutionT1S(1,:,5)', SolutionT1S(2,:,5)', SolutionT1S(3,:,5)', SolutionT1S(4,:,5)', SolutionT1S(5,:,5)', SolutionT1S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.7','0.965','1.3','1.7','2.0','2.5'}; set(gca,'XTickLabel',Labels); xlabel('T_{1S} (s)'); hline(0.2, 'm--');
% subplot(2,3,3); boxplot([Solution(:,5), SolutionT2F(1,:,5)', SolutionT2F(2,:,5)', SolutionT2F(3,:,5)', SolutionT2F(4,:,5)', SolutionT2F(5,:,5)', SolutionT2F(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','2','12','20','30','40','45'}; set(gca,'XTickLabel',Labels); xlabel('T_{2F} (ms)'); hline(0.2, 'm--');
% subplot(2,3,4); boxplot([Solution(:,5), SolutionT2S(1,:,5)', SolutionT2S(2,:,5)', SolutionT2S(3,:,5)', SolutionT2S(4,:,5)', SolutionT2S(5,:,5)', SolutionT2S(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(19); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','75','90','125','150','175','200'}; set(gca,'XTickLabel',Labels); xlabel('T_{2S} (ms)'); hline(0.2, 'm--');
% subplot(2,3,5); boxplot([Solution(:,5), SolutionkFS(1,:,5)', SolutionkFS(2,:,5)', SolutionkFS(3,:,5)', SolutionkFS(4,:,5)', SolutionkFS(5,:,5)', SolutionkFS(6,:,5)']); 
% a = get(get(gca,'children'),'children'); box = a(18); set(box,'Color','k')
% ylabel('MWF'); Labels = {'All','0.5','5','8','10','15','20'}; set(gca,'XTickLabel',Labels); xlabel('k_{FS} (s^{-1})'); hline(0.2, 'm--');

% figure(3);
% subplot(2,1,1); boxplot([Solution_FixedT2F(1,:,5)',Solution_FixedT2F(2,:,5)',Solution_FixedT2F(3,:,5)',Solution_FixedT2F(4,:,5)',Solution_FixedT2F(5,:,5)',Solution_FixedT2F(6,:,5)',Solution_FixedT2F(7,:,5)']);
% xlabel('SNR'), ylabel('MWF'); Labels = {'20','30','50','80','100','150','200'}; set(gca,'XTickLabel',Labels); hline(0.2, 'm--'); title('Fixed T_{2F}');
% subplot(2,1,2); boxplot([Solution_FixedT2S(1,:,5)',Solution_FixedT2S(2,:,5)',Solution_FixedT2S(3,:,5)',Solution_FixedT2S(4,:,5)',Solution_FixedT2S(5,:,5)',Solution_FixedT2S(6,:,5)',Solution_FixedT2S(7,:,5)']);
% xlabel('SNR'), ylabel('MWF'); Labels = {'20','30','50','80','100','150','200'}; set(gca,'XTickLabel',Labels); hline(0.2, 'm--'); title('Fixed T_{2S}');

figure(4);
boxplot([Solution_DB(1,:,5)', Solution_DD(1,:,5)', Solution_DW(1,:,5)', Solution_DZ(1,:,5)',Solution_RB(1,:,5)', Solution_RD(1,:,5)', Solution_RW(1,:,5)', Solution_RZ(1,:,5)']);
ylabel('MWF'); Labels = {'D+B','D+D','D+W','D+Z','R+B','R+D','R+W','R+Z'}; set(gca,'XTickLabel',Labels); hline(0.2, 'm--');