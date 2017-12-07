%% Fitting a mcDESPOT model to 2-pool MT data.

clear all

Runs = 25;
options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8);

% Define sequence and ground-truth tissue parameters.
TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_W = 1.15; T2_W = 0.08; M0_B = 0.25; k_WB = 9; T1_B = 0.4;

T1S_app = zeros(Runs,1); T2S_app = zeros(Runs,1); T1F_app = zeros(Runs,1); T2F_app = zeros(Runs,1); M0F_app = zeros(Runs,1); M0S_app = zeros(Runs,1); kFS_app = zeros(Runs,1); kSF_app = zeros(Runs,1);
T1S_rand = zeros(Runs,1); T2S_rand = zeros(Runs,1);
T1F_rand = zeros(Runs,1); T2F_rand = zeros(Runs,1);
M0F_rand = zeros(Runs,1); kFS_rand = zeros(Runs,1);
M0S_rand = zeros(Runs,1); kSF_rand = zeros(Runs,1);

Fval = zeros(Runs,1); Exitflag = zeros(Runs,1); Grad = zeros(8, Runs); Hessian = zeros(8, 8, Runs);

% GT-signals. Change between CS and US in functions.
SPGR_Data = TwoPoolMT_SPGR(FA_SPGR,TR_SPGR,'T1_W',T1_W,'T1_B',T1_B,'M0_B',M0_B,'k_WB',k_WB);
SSFP_Data = TwoPoolMT_SSFP(FA_SSFP,TR_SSFP,'T1_W',T1_W,'T1_B',T1_B,'T2_W',T2_W,'M0_B',M0_B,'k_WB',k_WB);
SSFP_Data_180 = TwoPoolMT_SSFP_180(FA_SSFP,TR_SSFP,'T1_W',T1_W,'T1_B',T1_B,'T2_W',T2_W,'M0_B',M0_B,'k_WB',k_WB);
Data = [SPGR_Data ; SSFP_Data ; SSFP_Data_180];

% T1S_LB = 0.2; T1S_UB = 1.25; T2S_LB = 0.04; T2S_UB = 0.15;
% T1F_LB = 0.1; T1F_UB = 0.5; T2F_LB = 0.01; T2F_UB = 0.03;
% M0F_LB = 0.05; M0F_UB = 0.3; M0S_LB = 0.05; M0S_UB = 0.8;
% kFS_LB = 1; kFS_UB = 20; kSF_LB = 1; kSF_UB = 10; % Original bounds.

T1S_LB = 0; T1S_UB = 3; T2S_LB = 0; T2S_UB = 0.5;
T1F_LB = 0; T1F_UB = 1; T2F_LB = 0; T2F_UB = 0.3;
M0F_LB = 0; M0F_UB = 0.5; M0S_LB = 0; M0S_UB = 0.8;
kFS_LB = 0; kFS_UB = 40; kSF_LB = 0; kSF_UB = 20;

tic
for ii = 1:Runs
    
    disp(ii)
    
    % Randomly specify parameter values for fitting starting point.
    T1S_rand(ii) = T1S_LB + (T1S_UB - T1S_LB) .* rand(1,1); T2S_rand(ii) = T2S_LB + (T2S_UB - T2S_LB) .* rand(1,1);
    T1F_rand(ii) = T1F_LB + (T1F_UB - T1F_LB) .* rand(1,1); T2F_rand(ii) = T2F_LB + (T2F_UB - T2F_LB) .* rand(1,1);
    M0F_rand(ii) = M0F_LB + (M0F_UB - M0F_LB) .* rand(1,1); kFS_rand(ii) = kFS_LB + (kFS_UB - kFS_LB) .* rand(1,1);
    M0S_rand(ii) = M0S_LB + (M0S_UB - M0S_LB) .* rand(1,1); kSF_rand(ii) = kSF_LB + (kSF_UB - kSF_LB) .* rand(1,1);
    x0 = [T1S_rand(ii) T1F_rand(ii) M0F_rand(ii) M0S_rand(ii) kFS_rand(ii) kSF_rand(ii) T2S_rand(ii) T2F_rand(ii)];
    
    Sig_SPGR = @(x)(SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',x(1),'T1_F',x(2),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
    Sig_SSFP = @(x)(SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
    Sig_SSFP180 = @(x)(SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',x(1),'T2_S',x(7),'T1_F',x(2),'T2_F',x(8),'M0_F',x(3),'M0_S',x(4),'k_FS',x(5),'k_SF',x(6)));
    Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
    CF = @(x)norm(Sig(x) - Data)^2;
    [Sol,Fval(ii),Exitflag(ii),Output,Lambda,Grad(:,ii),Hessian(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1S_LB T1F_LB M0F_LB M0S_LB kFS_LB kSF_LB T2S_LB T2F_LB],[T1S_UB T1F_UB M0F_UB M0S_UB kFS_UB kSF_UB T2S_UB T2F_UB], [], options);
    
    % Remove non-convergent solutions.
    if Exitflag(ii) == 1 %~= 0
        
        T1S_app(ii) = Sol(1); T1F_app(ii) = Sol(2); M0F_app(ii) = Sol(3); M0S_app(ii) = Sol(4); kFS_app(ii) = Sol(5); kSF_app(ii) = Sol(6); T2S_app(ii) = Sol(7); T2F_app(ii) = Sol(8);
        
    end
    
end
toc

[~,Index] = min(Fval);

T1S_final = T1S_app(Index); T1F_final = T1F_app(Index); M0F_final = M0F_app(Index); M0S_final = M0S_app(Index); kFS_final = kFS_app(Index); kSF_final = kSF_app(Index); T2S_final = T2S_app(Index); T2F_final = T2F_app(Index);

% Signal plots. Change title and subplot indices for CS or US.
figure(1)
subplot(3,2,1)
plot(rad2deg(FA_SPGR), (SPGR_steady_state_M0(FA_SPGR,TR_SPGR,'T1_S',T1S_final,'T1_F',T1F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SPGR), SPGR_Data, 'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('SPGR Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
title(['US (T_{1F}^{app} = ' num2str(round(T1F_final,2)) 's,' ' T_{2F}^{app} = ' num2str(round(T2F_final,2,'significant')) 's,' ' T_{1S}^{app} = ' num2str(round(T1S_final,2)) 's,' ' T_{2S}^{app} = ' num2str(round(T2S_final,2)) 's,' ' M_{0F}^{app} = ' num2str(round(M0F_final,2)) ',' ' M_{0S}^{app} = ' num2str(round(M0S_final,2)) ',' ' k_{FS}^{app} = ' num2str(round(kFS_final)) 's^{-1},' ' k_{SF}^{app} = ' num2str(round(kSF_final)) 's^{-1}' ')']);
ax = gca;
ax.TitleFontSizeMultiplier = 0.6;

subplot(3,2,5)
plot(rad2deg(FA_SSFP), (SSFP_steady_state_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{0} Signal','Fontsize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)
ll = legend('Fitted mcDESPOT Signal', 'Two-Pool MT Data');
ll.FontSize = 12;
set(ll.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.5;.5;.5;.1]));

subplot(3,2,3)
plot(rad2deg(FA_SSFP), (SSFP_steady_state_180_M0(FA_SSFP,TR_SSFP,'T1_S',T1S_final,'T2_S',T2S_final,'T1_F',T1F_final,'T2_F',T2F_final,'M0_F',M0F_final,'M0_S',M0S_final,'k_FS',kFS_final,'k_SF',kSF_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data_180,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',12); ylabel('bSSFP_{180} Signal','FontSize',12)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 12)

%% Fitting a single-pool to 3-pool data.

clear all

Runs = 10;
options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'OptimalityTolerance', 5e-8);

% Define sequence and ground-truth tissue parameters.
TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3;
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T1_F = 0.4; T1_B = 1; T2_S = 0.08; T2_F = 0.02; M0_F = 0.25; M0_S = 0.55; M0_B = 0.2; k_FS = 9; k_SB = 5; k_FB = 5;

% GT-signals. Change between CS and US here.
SPGR_Data = CS_SPGR_steady_state_MT_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data = CS_SSFP_steady_state_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = CS_SSFP_steady_state_180_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP_Data; SSFP_Data_180];

T1_LB = 0.1; T1_UB = 1.25; T2_LB = 0.01; T2_UB = 0.15; M0_LB = 0.05; M0_UB = 1.5;

T1_app = zeros(Runs,1); T2_app = zeros(Runs,1); M0_app = zeros(Runs,1);
T1_rand = zeros(Runs,1); T2_rand = zeros(Runs,1); M0_rand = zeros(Runs,1);

tic
for ii = 1:Runs
    
    disp(ii)
    
    % Randomly specify parameter values for fitting starting point.
    T1_rand(ii) = T1_LB + (T1_UB - T1_LB) .* rand(1,1); T2_rand(ii) = T2_LB + (T2_UB - T2_LB) .* rand(1,1);
    M0_rand(ii) = M0_LB + (M0_UB - M0_LB) .* rand(1,1);
    x0 = [T1_rand(ii) M0_rand(ii) T2_rand(ii)];
    
    Sig_SPGR = @(x)(SPGR_SP(FA_SPGR,TR_SPGR,'T1',x(1),'M0',x(2)));
    Sig_SSFP = @(x)(SSFP_SP(FA_SSFP,TR_SSFP,'T1',x(1),'T2',x(3),'M0',x(2)));
    Sig_SSFP180 = @(x)(SSFP_SP_180(FA_SSFP,TR_SSFP,'T1',x(1),'T2',x(3),'M0',x(2)));
    Sig = @(x)([Sig_SPGR(x) ; Sig_SSFP(x) ; Sig_SSFP180(x)]);
    CF = @(x)norm(Sig(x) - Data).^2;
    [Sol,Fval(ii),Exitflag(ii),Output,Lambda,Grad(:,ii),Hessian(:,:,ii)] = fmincon(CF,x0,[],[],[],[],[T1_LB M0_LB T2_LB],[T1_UB M0_UB T2_UB], [], options);
    
    % Remove non-convergent solutions.
    if Exitflag(ii) ~= 0
        
        T1_app(ii) = Sol(1); M0_app(ii) = Sol(2); T2_app(ii) = Sol(3);
    
    end
    
end
toc
[~,Index] = min(Fval);
T1_final = T1_app(Index); M0_final = M0_app(Index); T2_final = T2_app(Index);

% Signal plots. Change title and subplot indices for CS or US.
figure(1)
hold on
subplot(3,2,2)
plot(rad2deg(FA_SPGR), (SPGR_SP(FA_SPGR,TR_SPGR,'T1',T1_final,'M0',M0_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SPGR), SPGR_Data, 'ro','LineWidth',2)
tt = title(['CS (T_1^{app} = ' num2str(round(T1_final,2)) 's,' ' T_2^{app} = ' num2str(round(T2_final,2)) 's,' ' M_0^{app} = ' num2str(round(M0_final,2)) ')'],'FontWeight','normal');
xlabel('FA [deg]','FontSize',18); ylabel('SPGR Signal','Fontsize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
subplot(3,2,6)
plot(rad2deg(FA_SSFP), (SSFP_SP(FA_SSFP,TR_SSFP,'T1',T1_final,'T2',T2_final,'M0',M0_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{0} Signal','Fontsize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
ll = legend('Fitted Single-Pool Signal', 'Three-Pool Data');
ll.FontSize = 18;
set(ll.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.5;.5;.5;.1]));
subplot(3,2,4)
plot(rad2deg(FA_SSFP), (SSFP_SP_180(FA_SSFP,TR_SSFP,'T1',T1_final,'T2',T2_final,'M0',M0_final)),'k--','LineWidth',2)
hold on
plot(rad2deg(FA_SSFP), SSFP_Data_180,'ro','LineWidth',2)
xlabel('FA [deg]','FontSize',18); ylabel('bSSFP_{180} Signal','FontSize',18)
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
