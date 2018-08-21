%% mcDESPOT processing using simulated data.
close all
clear all

% Specify options for fmincon.
%rng('default')
options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8);

% Specify number of repeats.
Runs = 10; dMzB = 0; Sigma = 1/500;

% Define sequence and ground-truth tissue parameters.
TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; % [4.8 and 5.6]
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T2_S = 0.08; T1_F = 0.4; T2_F = 0.02; M0_F = 0.25; M0_S = 0.55; k_FS = 9; T1_B = 1; T2_B = 8e-6; M0_B = 0.2; k_SB = 5; k_FB = 5; % [Liu, T1s = 1.15, fb = 0.2]

% Initialise loop vectors.
T1S_app = zeros(Runs,length(M0_B)); T1F_app = zeros(Runs,length(M0_B)); M0F_app = zeros(Runs,length(M0_B)); M0S_app = zeros(Runs,length(M0_B)); kFS_app = zeros(Runs,length(M0_B)); T2S_app = zeros(Runs,length(M0_B)); T2F_app = zeros(Runs,length(M0_B)); kSF_app = zeros(Runs,length(M0_B));
T1S_final = zeros(length(M0_B),1); T1F_final = zeros(length(M0_B),1); T2S_final = zeros(length(M0_B),1); T2F_final = zeros(length(M0_B),1); M0F_final = zeros(length(M0_B),1); M0S_final = zeros(length(M0_B),1); kFS_final = zeros(length(M0_B),1); kSF_final = zeros(length(M0_B),1);
Error_T1S = zeros(length(M0_B),1); Error_T2S = zeros(length(M0_B),1);
Error_T1F = zeros(length(M0_B),1); Error_T2F = zeros(length(M0_B),1);
Error_M0F = zeros(length(M0_B),1); Error_kFS = zeros(length(M0_B),1);
Precision_T1S = zeros(length(M0_B),1); Precision_T2S = zeros(length(M0_B),1);
Precision_T1F = zeros(length(M0_B),1); Precision_T2F = zeros(length(M0_B),1);
Precision_M0F = zeros(length(M0_B),1); Precision_kFS = zeros(length(M0_B),1);

kSF_Prediction = zeros(length(M0_B),1); kFS_Prediction = zeros(length(M0_B),1);
R1F_Prediction = zeros(length(M0_B),1); R1S_Prediction = zeros(length(M0_B),1);
M0S_Prediction = zeros(length(M0_B),1); M0F_Prediction = zeros(length(M0_B),1);
T1S_Prediction = zeros(length(M0_B),1); T1F_Prediction = zeros(length(M0_B),1);
kSF_Prediction2 = zeros(length(M0_B),1); kFS_Prediction2 = zeros(length(M0_B),1);
R1F_Prediction2 = zeros(length(M0_B),1); R1S_Prediction2 = zeros(length(M0_B),1);
M0S_Prediction2 = zeros(length(M0_B),1); M0F_Prediction2 = zeros(length(M0_B),1);
T1S_Prediction2 = zeros(length(M0_B),1); T1F_Prediction2 = zeros(length(M0_B),1);

T1S_rand = zeros(Runs,length(M0_B)); T2S_rand = zeros(Runs,length(M0_B));
T1F_rand = zeros(Runs,length(M0_B)); T2F_rand = zeros(Runs,length(M0_B));
M0F_rand = zeros(Runs,length(M0_B)); kFS_rand = zeros(Runs,length(M0_B));
M0S_rand = zeros(Runs,length(M0_B)); kSF_rand = zeros(Runs,length(M0_B));
MWF = zeros(length(M0_B),1);

for ss = 1:length(M0_B)
    
    Fval = zeros(Runs,1); Exitflag = zeros(Runs,1); Grad = zeros(8, Runs); Hessian = zeros(8, 8, Runs);
    
    MWF(ss) = M0_F/(1-M0_B(ss));
    
    % Calulcate parameter predictions from derived expressions.
    R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S;
    G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
    T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SPGR);
    % From symbolic maths.
    kSF_Prediction = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    kFS_Prediction = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    R1F_Prediction = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    R1S_Prediction = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
    M0F_Prediction = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
    M0S_Prediction = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
    
    T1F_Prediction(ss) = 1/R1F_Prediction(ss); T1S_Prediction(ss) = 1/R1S_Prediction(ss);
    
    % From algebra.
    k_BS = k_SB*M0_S/M0_B; k_BF = k_FB*M0_F/M0_B; k_SF = k_FS*M0_F/M0_S;
    kFS_Prediction2(ss) = k_FS + ((k_BS*k_FB)/(R1_B + k_BS + k_BF + W));
    kSF_Prediction2(ss) = k_SF + ((k_BF*k_SB)/(R1_B + k_BS + k_BF + W));
    R1F_Prediction2(ss) = R1_F + k_FB - ((k_BS*k_FB + k_BF*k_FB)/(R1_B + k_BS + k_BF + W));
    R1S_Prediction2(ss) = R1_S + k_SB - ((k_BF*k_SB + k_BS*k_SB)/(R1_B + k_BS + k_BF + W));
    M0F_Prediction2(ss) = (M0_F*R1_F + (((M0_B*k_BF*R1_B)-dMzB*k_BF)/(R1_B + k_BS + k_BF + W)))/R1F_Prediction2(ss);
    M0S_Prediction2(ss) = (M0_S*R1_S + (((M0_B*k_BS*R1_B)-dMzB*k_BS)/(R1_B + k_BS + k_BF + W)))/R1S_Prediction2(ss);
    
    T1F_Prediction2(ss) = 1/R1F_Prediction2(ss); T1S_Prediction2(ss) = 1/R1S_Prediction2(ss);
    
    % GT-signals. dMzB = 0 for steady-state sequences.
    SPGR_Data = CS_SPGR_steady_state_MT_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    SSFP_Data = CS_SSFP_steady_state_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    SSFP_Data_180 = CS_SSFP_steady_state_180_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B(ss),'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
    Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];
    
    % Add noise to GT-signals.
    % M0 = 1; SNR = 100; NF = (0.1 * M0)/SNR;
    % SPGR_Data_Noisy = SPGR_Data + complex(NF * randn(size(SPGR_Data)), NF * randn(size(SPGR_Data)));
    % SSFP_Data_Noisy = SSFP_Data + complex(NF * randn(size(SSFP_Data)), NF * randn(size(SSFP_Data)));
    % SSFP_Data_Noisy_180 = SSFP_Data_180 + complex(NF * randn(size(SSFP_Data_180)), NF * randn(size(SSFP_Data_180)));
    % Concatenate signals for mcDESPOT processing. Choose noisy or noiseless.
    % Data = [SPGR_Data_Noisy ; SSFP_Data_Noisy ; SSFP_Data_Noisy_180];

    for ii = 1:Runs
        
        % Bounds of search-space.
        %T1S_LB = 0.2; T1S_UB = 1.5; T2S_LB = 0.04; T2S_UB = 0.15; % [T1s_LB = 0.3]
        %T1F_LB = 0.1; T1F_UB = 0.8; T2F_LB = 0.01; T2F_UB = 0.03;
        %M0F_LB = 0.001; M0F_UB = 0.35; M0S_LB = 0.1; M0S_UB = 0.5; % [M0S_LB = 0.3, M0S_UB = 0.8]
        %kFS_LB = (1/0.6); kFS_UB = 40; kSF_LB = (1/0.6); kSF_UB = 40; % WOOD BUT ALTERED LBs of T1s.
        
        T1S_LB = 0.2; T1S_UB = 1.25; T2S_LB = 0.04; T2S_UB = 0.15;
        T1F_LB = 0.1; T1F_UB = 0.5; T2F_LB = 0.01; T2F_UB = 0.03;
        M0F_LB = 0.05; M0F_UB = 0.3; M0S_LB = 0.05; M0S_UB = 0.8;
        kFS_LB = 1; kFS_UB = 20; kSF_LB = 1; kSF_UB = 10; % Tightened bounds.
        
        % Randomly specify parameter values for fitting starting point.
        T1S_rand(ii,ss) = T1S_LB + (T1S_UB - T1S_LB) .* rand(1,1); T2S_rand(ii,ss) = T2S_LB + (T2S_UB - T2S_LB) .* rand(1,1);
        T1F_rand(ii,ss) = T1F_LB + (T1F_UB - T1F_LB) .* rand(1,1); T2F_rand(ii,ss) = T2F_LB + (T2F_UB - T2F_LB) .* rand(1,1);
        M0F_rand(ii,ss) = M0F_LB + (M0F_UB - M0F_LB) .* rand(1,1); kFS_rand(ii,ss) = kFS_LB + (kFS_UB - kFS_LB) .* rand(1,1);
        M0S_rand(ii,ss) = M0S_LB + (M0S_UB - M0S_LB) .* rand(1,1); kSF_rand(ii,ss) = kSF_LB + (kSF_UB - kSF_LB) .* rand(1,1);
        x0 = [T1S_rand(ii,ss) T1F_rand(ii,ss) M0F_rand(ii,ss) M0S_rand(ii,ss) kFS_rand(ii,ss) kSF_rand(ii,ss) T2S_rand(ii,ss) T2F_rand(ii,ss)];
        %x0 = [T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F];
        %x0 = [T1_S T1_F M0_F M0_S k_FS k_SF T2_S T2_F];
        
        % Perform fitting of two-pool models to three-pool signals.
        disp(ii)
        % Non-MT.
        Sig_SPGR = @(x)(SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig_SSFP = @(x)(SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig_SSFP180 = @(x)(SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
        Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
        CF = @(x) norm(Sig(x) - Data)^2; %(1/(2*Sigma^2))
        [Sol,Fval(ii),Exitflag(ii),Output,Lambda,Grad(:,ii),Hessian(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1S_LB T1F_LB M0F_LB M0S_LB kFS_LB kSF_LB T2S_LB T2F_LB],[T1S_UB T1F_UB M0F_UB M0S_UB kFS_UB kSF_UB T2S_UB T2F_UB], [], options);
        
        % Remove non-convergent solutions.
        if Exitflag(ii) == 1
            
            T1S_app(ii,ss) = Sol(1); T1F_app(ii,ss) = Sol(2); M0F_app(ii,ss) = Sol(3); M0S_app(ii,ss) = Sol(4); kFS_app(ii,ss) = Sol(5); kSF_app(ii,ss) = Sol(6); T2S_app(ii,ss) = Sol(7); T2F_app(ii,ss) = Sol(8);
            
        end
        
    end
    
    % Average parameter values calculated across 'Runs'.
    %T1S_final(ss) = mean(nonzeros(T1S_app(:,ss))); T1F_final(ss) = mean(nonzeros(T1F_app(:,ss))); T2S_final(ss) = mean(nonzeros(T2S_app(:,ss))); T2F_final(ss) = mean(nonzeros(T2F_app(:,ss))); M0F_final(ss) = mean(nonzeros(M0F_app(:,ss))); M0S_final(ss) = mean(nonzeros(M0S_app(:,ss))); kFS_final(ss) = mean(nonzeros(kFS_app(:,ss))); kSF_final(ss) = mean(nonzeros(kSF_app(:,ss)));
    [MinFval,Index] = min(Fval);
    T1S_final(ss) = T1S_app(Index,ss); T1F_final(ss) = T1F_app(Index,ss); M0F_final(ss) = M0F_app(Index,ss); M0S_final(ss) = M0S_app(Index,ss); kFS_final(ss) = kFS_app(Index,ss); kSF_final(ss) = kSF_app(Index,ss); T2S_final(ss) = T2S_app(Index,ss); T2F_final(ss) = T2F_app(Index,ss);
    
    % Calculate accuracy and precision of parameter estimates from fmincon. MWF not used?
    Error_T1S(ss) = abs(((T1S_final(ss) - T1_S)/T1_S)*100); Error_T2S(ss) = abs(((T2S_final(ss) - T2_S)/T2_S)*100);
    Error_T1F(ss) = abs(((T1F_final(ss) - T1_F)/T1_F)*100); Error_T2F(ss) = abs(((T2F_final(ss) - T2_F)/T2_F)*100);
    Error_M0F(ss) = abs(((M0F_final(ss) - M0_F)/M0_F)*100); Error_kFS(ss) = abs(((kFS_final(ss) - k_FS)/k_FS)*100);
    Precision_T1S(ss) = std(T1S_app(:,ss))/(T1S_final(ss)); Precision_T2S(ss) = std(T2S_app(:,ss))/(T2S_final(ss));
    Precision_T1F(ss) = std(T1F_app(:,ss))/(T1F_final(ss)); Precision_T2F(ss) = std(T2F_app(:,ss))/(T2F_final(ss));
    Precision_M0F(ss) = std(M0F_app(:,ss))/(M0F_final(ss)); Precision_kFS(ss) = std(kFS_app(:,ss))/(kFS_final(ss));
    
    % Signal plots and histograms. Change functions if three-pool or noisy fit is used.
    figure(1)
    hold on
    subplot(3,2,1)
    plot(rad2deg(FA_SPGR), (SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_final(ss),'T1_F',T1F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k','LineWidth',2)
    hold on
    plot(rad2deg(FA_SPGR), SPGR_Data, 'bo','LineWidth',2)
    %plot(rad2deg(FA_SPGR), abs(SPGR_Data_Noisy), 'bo','LineWidth',2)
    plot(rad2deg(FA_SPGR),SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_Prediction(ss),'T1_F',T1F_Prediction(ss),'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth', 2)
    %plot(rad2deg(FA_SPGR),SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_Prediction2(ss),'T1_F',T1F_Prediction2(ss),'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth', 2)
    ll = legend('Fitted Signal', 'Three-Pool CS Data', 'Two-Pool Apparent Signal');
    ll.FontSize = 18;
    xlabel('FA [deg]','FontSize',18); ylabel('SPGR Signal','Fontsize',18)
    get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
    title('(b) Uncontrolled Saturation','FontSize',18)
    subplot(3,2,5)
    plot(rad2deg(FA_SSFP), (SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final(ss),'T2_S',T2S_final(ss),'T1_F',T1F_final(ss),'T2_F',T2F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k','LineWidth',2)
    hold on
    plot(rad2deg(FA_SSFP), SSFP_Data,'bo','LineWidth',2)
    %plot(rad2deg(FA_SSFP), abs(SSFP_Data_Noisy),'bo','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_M0(FA_SSFP, TR_SSFP, 'T1_S',T1S_Prediction(ss),'T2_S',T2_S,'T1_f',T1F_Prediction(ss),'T2_F',T2_F,'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth',2)
    %plot(rad2deg(FA_SSFP), SSFP_steady_state_M0(FA_SSFP, TR_SSFP, 'T1_S',T1S_Prediction2(ss),'T2_S',T2_S,'T1_f',T1F_Prediction2(ss),'T2_F',T2_F,'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth',2)
    ll = legend('Fitted Signal', 'Three-Pool CS Data','Two-Pool Apparent Signal');
    ll.FontSize = 18;
    xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{0} Signal','Fontsize',18)
    get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
    subplot(3,2,3)
    plot(rad2deg(FA_SSFP), (SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final(ss),'T2_S',T2S_final(ss),'T1_F',T1F_final(ss),'T2_F',T2F_final(ss),'M0_F',M0F_final(ss),'M0_S',M0S_final(ss),'k_FS',kFS_final(ss),'k_SF',kSF_final(ss))),'k','LineWidth',2)
    hold on
    plot(rad2deg(FA_SSFP), SSFP_Data_180,'bo','LineWidth',2)
    %plot(rad2deg(FA_SSFP), abs(SSFP_Data_Noisy_180),'bo','LineWidth',2)
    plot(rad2deg(FA_SSFP), SSFP_steady_state_180_M0(FA_SSFP, TR_SSFP, 'T1_S',T1S_Prediction(ss),'T2_S',T2_S,'T1_F',T1F_Prediction(ss),'T2_F',T2_F,'M0_F',M0F_Prediction(ss),'M0_S',M0S_Prediction(ss),'k_FS',kFS_Prediction(ss),'k_SF',kSF_Prediction(ss)), 'r--','LineWidth',2)
    %plot(rad2deg(FA_SSFP), SSFP_steady_state_180_M0(FA_SSFP, TR_SSFP, 'T1_S',T1S_Prediction2(ss),'T2_S',T2_S,'T1_F',T1F_Prediction2(ss),'T2_F',T2_F,'M0_F',M0F_Prediction2(ss),'M0_S',M0S_Prediction2(ss),'k_FS',kFS_Prediction2(ss),'k_SF',kSF_Prediction2(ss)), 'mo','LineWidth',2)
    ll = legend('Fitted Signal', 'Three-Pool CS Data','Two-Pool Apparent Signal');
    ll.FontSize = 18;
    xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{180} Signal','FontSize',18)
    get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
    
    %     figure(10)
    %     subplot(2,4,1); His3 = histogram(nonzeros(T1S_app(:,ss)));
    %     hold on; vline(T1_S,'k','GT')
    %     xlabel('T_{1S}'); ylabel('Count');
    %     His3.FaceColor = 'r'; His3.EdgeColor = 'k';
    %     subplot(2,4,2); His4 = histogram(nonzeros(T1F_app(:,ss)));
    %     hold on; vline(T1_F,'k','GT')
    %     xlabel('T_{1F}'); ylabel('Count');
    %     His4.FaceColor = 'r'; His4.EdgeColor = 'k';
    %     subplot(2,4,3); His5 = histogram(nonzeros(M0F_app(:,ss)));
    %     hold on; vline(M0_F,'k','GT')
    %     xlabel('M0_{F}'); ylabel('Count');
    %     His5.FaceColor = 'r'; His5.EdgeColor = 'k';
    %     subplot(2,4,4); His6 = histogram(nonzeros(kFS_app(:,ss)));
    %     hold on; vline(k_FS,'k','GT')
    %     xlabel('k_{FS}'); ylabel('Count');
    %     His6.FaceColor = 'r'; His6.EdgeColor = 'k';
    %     subplot(2,4,5); His7 = histogram(nonzeros(T2S_app(:,ss)));
    %     hold on; vline(T2_S,'k','GT')
    %     xlabel('T_{2S}'); ylabel('Count');
    %     His7.FaceColor = 'r'; His7.EdgeColor = 'k';
    %     subplot(2,4,6); His8 = histogram(nonzeros(T2F_app(:,ss)));
    %     hold on; vline(T2_F,'k','GT')
    %     xlabel('T_{2F}'); ylabel('Count');
    %     His8.FaceColor = 'r'; His8.EdgeColor = 'k';
    %     subplot(2,4,7); His9 = histogram(nonzeros(M0S_app(:,ss)));
    %     hold on; vline(M0_S,'k','GT')
    %     xlabel('M_{0S}'); ylabel('Count');
    %     His9.FaceColor = 'r'; His9.EdgeColor = 'k';
    %     subplot(2,4,8); His10 = histogram(nonzeros(kSF_app(:,ss)));
    %     hold on; vline(k_SF,'k','GT')
    %     xlabel('k_{SF}'); ylabel('Count');
    %     His10.FaceColor = 'r'; His10.EdgeColor = 'k';
    
end