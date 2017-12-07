%% mcDESPOT processing using simulated data.

%close all
%clear all

% Specify options for fmincon.
%rng('default')
options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, 'OptimalityTolerance', 5e-8);

% Specify number of repeats.
Runs = 1000; dMzB = 0; Sigma = 1/5000; %linspace((1/5000),(1/200),Runs);

% Define sequence and ground-truth tissue parameters.
TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; % [4.8 and 5.6]
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T2_S = 0.08; T1_F = 0.4; T2_F = 0.02; M0_F = 0.25; k_FS = 9;

% Initialise loop vectors.
T1S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); M0F_app = zeros(Runs,1); M0S_app = zeros(Runs,1); kFS_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T2F_app = zeros(Runs,1); kSF_app = zeros(Runs,1);
T1S_rand = zeros(Runs,1); T2S_rand = zeros(Runs,1);
T1F_rand = zeros(Runs,1); T2F_rand = zeros(Runs,1);
M0F_rand = zeros(Runs,1); kFS_rand = zeros(Runs,1);
M0S_rand = zeros(Runs,1); kSF_rand = zeros(Runs,1);

Fval = zeros(Runs,1); Exitflag = zeros(Runs,1); Grad = zeros(6, Runs); Hessian = zeros(6, 6, Runs);

% GT-signals. dMzB = 0 for steady-state sequences.
SPGR_Data = SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data = SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
SSFP_Data_180 = SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'M0_F',M0_F,'k_FS',k_FS);
Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];

T1S_LB = 0.9; T1S_UB = 1.5; T2S_LB = 0.04; T2S_UB = 0.15;
T1F_LB = 0.3; T1F_UB = 0.8; T2F_LB = 0.01; T2F_UB = 0.03;
M0F_LB = 0.001; M0F_UB = 0.35; kFS_LB = (1/0.6); kFS_UB = 40;
%x0 = [T1_S T1_F M0_F k_FS T2_S T2_F];

Eigenvalues = zeros(6,1,Runs); Ratio = zeros(Runs, 1);

disp (['Estimated Run Time: ' num2str(1.352192*Runs) ' seconds']);
tic
for ii = 1:Runs
    
    % Randomly specify parameter values for fitting starting point.
    T1S_rand(ii) = T1S_LB + (T1S_UB - T1S_LB) .* rand(1,1); T2S_rand(ii) = T2S_LB + (T2S_UB - T2S_LB) .* rand(1,1);
    T1F_rand(ii) = T1F_LB + (T1F_UB - T1F_LB) .* rand(1,1); T2F_rand(ii) = T2F_LB + (T2F_UB - T2F_LB) .* rand(1,1);
    M0F_rand(ii) = M0F_LB + (M0F_UB - M0F_LB) .* rand(1,1); kFS_rand(ii) = kFS_LB + (kFS_UB - kFS_LB) .* rand(1,1);
    x0 = [T1S_rand(ii) T1F_rand(ii) M0F_rand(ii) kFS_rand(ii) T2S_rand(ii) T2F_rand(ii)];
    
    % Perform fitting of two-pool models to three-pool signals.
    disp(ii)
    % Non-MT.
    Sig_SPGR = @(x)(SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'k_FS',x(4)));
    Sig_SSFP = @(x)(SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
    Sig_SSFP180 = @(x)(SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(5),'T1_F',x(2),'T2_F',x(6),'M0_F',x(3),'k_FS',x(4)));
    Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
    CF = @(x)(1/(2*Sigma^2))*norm(Sig(x) - Data).^2; % Sigma(ii)
    [Sol,Fval(ii),Exitflag(ii),Output,Lambda,Grad(:,ii),Hessian(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1S_LB T1F_LB M0F_LB kFS_LB T2S_LB T2F_LB],[T1S_UB T1F_UB M0F_UB kFS_UB T2S_UB T2F_UB], [], options);
    
    % Remove non-convergent solutions.
    %if Exitflag(ii) == 1
        
        T1S_app(ii) = Sol(1); T1F_app(ii) = Sol(2); M0F_app(ii) = Sol(3); kFS_app(ii) = Sol(4); T2S_app(ii) = Sol(5); T2F_app(ii) = Sol(6);
        
        Eigenvalues(:,:,ii) = eig(Hessian(:,:,ii));
        Ratio(ii) = Eigenvalues(6,1,ii)/Eigenvalues(5,1,ii);
    
    %end
    
end
toc
[MinFval,Index] = min(Fval);
T1S_final = T1S_app(Index); T1F_final = T1F_app(Index); M0F_final = M0F_app(Index); kFS_final = kFS_app(Index); T2S_final = T2S_app(Index); T2F_final = T2F_app(Index);

% Signal plots and histograms. Change functions if three-pool or noisy fit is used.
figure(1)
hold on
subplot(3,2,2)
plot(rad2deg(FA_SPGR), (SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'k_FS',kFS_final)),'k','LineWidth',2)
hold on
plot(rad2deg(FA_SPGR), SPGR_Data, 'bo','LineWidth',2)
ll = legend('Fitted Signal', 'Three-Pool CS Data');
ll.FontSize = 18;
xlabel('FA [deg]','FontSize',18); ylabel('SPGR Signal','Fontsize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
title('(b) Controlled Saturation','FontSize',18)
subplot(3,2,6)
plot(rad2deg(FA_SSFP), (SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'k_FS',kFS_final)),'k','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data,'bo','LineWidth',2)
ll = legend('Fitted Signal', 'Three-Pool CS Data');
ll.FontSize = 18;
xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{0} Signal','Fontsize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
subplot(3,2,4)
plot(rad2deg(FA_SSFP), (SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'k_FS',kFS_final)),'k','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data_180,'bo','LineWidth',2)
ll = legend('Fitted Signal', 'Three-Pool CS Data');
ll.FontSize = 18;
xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{180} Signal','FontSize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)

% for ii = 1:1000
%     if Ratio(ii) > 1000
%         Ratio(ii) = 1000;
%     end
% end

figure(2); H = histogram(Ratio, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha',0.5); xlabel('Eigenvalue Ratio','FontSize',16); ylabel('Count','FontSize',18);
% hold on;
% histogram(Ratio(201:400), 'EdgeColor', 'k', 'FaceColor', 'g', 'FaceAlpha',0.5);
% histogram(Ratio(401:600), 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha',0.5);
% histogram(Ratio(601:800), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha',0.5);
% histogram(Ratio(801:1000), 'EdgeColor', 'k', 'FaceColor', 'y', 'FaceAlpha',0.5);
% l = legend();

disp (['Ratio Min: ' num2str(min(Ratio))]);
disp (['Ratio Max: ' num2str(max(Ratio))]);