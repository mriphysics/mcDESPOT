clear all

TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; % Make different?
FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); 
FA_SSFP = deg2rad([2 5 10 15 20 30 40 50]);
T1_S = 1.15; T1_F = 0.4; T1_B = 1;
T2_S = 0.08; T2_F = 0.02;
M0_F = 0.25; M0_S = 0.55; M0_B = 0.2;
k_FS = 9; k_SB = 5; k_FB = 5; k_SF = (M0_F*k_FS)/M0_S;
dMzB = 0;

% GT-signals.
SPGR_Data = CS_SPGR_steady_state_MT_M0(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data = CS_SSFP_steady_state_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = CS_SSFP_steady_state_180_MT_M0(FA_SSFP,TR_SSFP,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
Data = [SPGR_Data ; SSFP_Data; SSFP_Data_180];

% Calulcate parameter predictions from derived expressions.
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR_SPGR);

kSF_Prediction = (M0_F*k_FS)/M0_S + (M0_F*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
kFS_Prediction = k_FS + (M0_S*k_FB*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1F_Prediction = (M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
R1S_Prediction = (M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB)/(M0_B*W + M0_F*k_FB + M0_S*k_SB + M0_B*R1_B);
M0F_Prediction = (M0_F*(M0_B*R1_B*R1_F - dMzB*k_FB + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB))/(M0_B*R1_B*R1_F + M0_B*R1_F*W + M0_B*R1_B*k_FB + M0_F*R1_F*k_FB + M0_S*R1_F*k_SB + M0_B*W*k_FB);
M0S_Prediction = (M0_S*(M0_B*R1_B*R1_S - dMzB*k_SB + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB))/(M0_B*R1_B*R1_S + M0_B*R1_S*W + M0_F*R1_S*k_FB + M0_B*R1_B*k_SB + M0_S*R1_S*k_SB + M0_B*W*k_SB);
T1F_Prediction = 1/R1F_Prediction;
T1S_Prediction = 1/R1S_Prediction;

nSamples = 500;
Downsample = 10;
Runs = 4;
Sigma = 1/300; %SNR = 500.

disp (['Estimated Run Time: ' num2str((nSamples*Downsample/100)) ' seconds']);

% Upper = [1.5 0.75 0.5 0.75 20 10 0.2 0.1];
% Lower = [0 0 0 0 0 0 0 0];
Upper = [0.4 0.25 0.35 0.75 40 40 0.15 0.03];
Lower = [0.2 0.001 0.001 0.001 (1/0.6) (1/0.6) 0.04 0.01]; % MODIFIED.

Scale = 0.1*diag([0.001 0.001 0.0001 0.0001 0.1 0.1 0.001 0.001]);

Start = [T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F];

Proprng = @(x) mvnrnd(x,Scale,1);

lnP = @(x) logpdf(x, Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma);

tic

% Perform MHM.
Smpl = zeros(nSamples, 8, Runs); Starts = zeros(Runs,8); lnP_val = zeros(size(Smpl,1),1);
for ii = 1:Runs
    while 1 
      x = Proprng(Start);
      if isfinite(lnP(x)), Starts(ii,:) = x; break; end
    end
end

parfor ii = 1:Runs
    [Smpl(:,:,ii), Rate(ii)] = mhsample(Start, nSamples,'logpdf',lnP,'proprnd',Proprng,'symmetric',true,'thin', Downsample);
end
toc

Smpl = reshape(permute(Smpl, [1,3,2]),nSamples*Runs,8);
for ii=1:size(Smpl,1)
    %lnP_val(ii) = lnP (Smpl(ii,:));
    lnP_val(ii) = logpdf(Smpl(ii,:), Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma);
end

%% Plot results.

figure(1);
clf;
subplot(4,2,1); plot(Smpl(:,1)); title('T_{1S}')
hold on; hline(T1S_Prediction, 'k', 'Pred'); hline(T1_S, 'r', '3P')
subplot(4,2,2); plot(Smpl(:,2)); title('T_{1F}')
hold on; hline(T1F_Prediction, 'k', 'Pred'); hline(T1_F, 'r', '3P')
subplot(4,2,3); plot(Smpl(:,4)); title('M_{0S}')
hold on; hline(M0S_Prediction, 'k', 'Pred'); hline(M0_S, 'r', '3P')
subplot(4,2,4); plot(Smpl(:,3)); title('M_{0F}')
hold on; hline(M0F_Prediction, 'k', 'Pred'); hline(M0_F, 'r', '3P')
subplot(4,2,5); plot(Smpl(:,5)); title('k_{FS}')
hold on; hline(kFS_Prediction, 'k', 'Pred'); hline(k_FS, 'r', '3P')
subplot(4,2,6); plot(Smpl(:,6)); title('k_{SF}')
hold on; hline(kSF_Prediction, 'k', 'Pred'); hline(k_SF, 'r', '3P')
subplot(4,2,7); plot(Smpl(:,7)); title('T_{2S}')
hold on; hline(T2_S, 'k', '3P')
subplot(4,2,8); plot(Smpl(:,8)); title('T_{2F}')
hold on; hline(T2_F, 'k', '3P')

figure(2);
clf
subplot(4,2,1); h = histfit(Smpl(:,1)); h(1).FaceColor = [.8 .8 1]; title('T_{1S}')
hold on; vline(T1S_Prediction, 'k', 'Pred'); %vline(T1_S, 'r', '3P')
subplot(4,2,2); h = histfit(Smpl(:,2)); h(1).FaceColor = [.8 .8 1]; title('T_{1F}')
hold on; vline(T1F_Prediction, 'k', 'Pred'); %vline(T1_F, 'r', '3P')
subplot(4,2,3); h = histfit(Smpl(:,4)); h(1).FaceColor = [.8 .8 1]; title('M_{0S}')
hold on; vline(M0S_Prediction, 'k', 'Pred'); %vline(M0_S, 'r', '3P')
subplot(4,2,4); h = histfit(Smpl(:,3)); h(1).FaceColor = [.8 .8 1]; title('M_{0F}')
hold on; vline(M0F_Prediction, 'k', 'Pred'); %vline(M0_F, 'r', '3P')
subplot(4,2,5); h = histfit(Smpl(:,5)); h(1).FaceColor = [.8 .8 1]; title('k_{FS}')
hold on; vline(kFS_Prediction, 'k', 'Pred'); %vline(k_FS, 'r', '3P')
subplot(4,2,6); h = histfit(Smpl(:,6)); h(1).FaceColor = [.8 .8 1]; title('k_{SF}')
hold on; vline(kSF_Prediction, 'k', 'Pred'); %vline(k_SF, 'r', '3P')
subplot(4,2,7); h = histfit(Smpl(:,7)); h(1).FaceColor = [.8 .8 1]; title('T_{2S}')
hold on; vline(T2_S, 'k', '3P')
subplot(4,2,8); h = histfit(Smpl(:,8)); h(1).FaceColor = [.8 .8 1]; title('T_{2F}')
hold on; vline(T2_F, 'k', '3P')

figure(3); % Change depending on parameter space to visualise.
clf
h = plot3(Smpl(:,1), Smpl(:,2), Smpl(:,5), '.'); % T1S, T1F, kFS.
set(h, 'MarkerSize', 5)

%% Direct sampling of cost-function in 1D.

T1S_Vector = linspace(0.2,0.3,1000);
T1F_Vector = linspace(0.14,0.22,1000);
M0F_Vector = linspace(0.08,0.16,1000);
M0S_Vector = linspace(0.1,0.18,1000);
kFS_Vector = linspace(10,11,1000);
kSF_Vector = linspace(4.3,5.1,1000);
T2S_Vector = linspace(0.06,0.1,1000);
T2F_Vector = linspace(0.01,0.03,1000);

P_T1S = zeros(length(T1S_Vector),1);
P_T1F = zeros(length(T1S_Vector),1);
P_M0F = zeros(length(T1S_Vector),1);
P_M0S = zeros(length(T1S_Vector),1);
P_kFS = zeros(length(T1S_Vector),1);
P_kSF = zeros(length(T1S_Vector),1);
P_T2S = zeros(length(T1S_Vector),1);
P_T2F = zeros(length(T1S_Vector),1);

tic
for ii = 1:length(T1S_Vector)

    P_T1S(ii) = (logpdf([T1S_Vector(ii) T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_T1F(ii) = (logpdf([T1S_Prediction T1F_Vector(ii) M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_M0F(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_M0S(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Vector(ii) kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_kFS(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Vector(ii) kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_kSF(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Vector(ii) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_T2S(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2S_Vector(ii) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
    P_T2F(ii) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2F_Vector(ii)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));

end
toc

figure(4);
subplot(2,4,1); plot(T1S_Vector,exp(P_T1S),'k'); xlabel('T_{1S}'); ylabel('PDF Value'); hold on; vline(T1S_Prediction, 'r', 'Pred')
subplot(2,4,2); plot(T1F_Vector,exp(P_T1F),'k'); xlabel('T_{1F}'); ylabel('PDF Value'); hold on; vline(T1F_Prediction, 'r', 'Pred')
subplot(2,4,4); plot(M0F_Vector,exp(P_M0F),'k'); xlabel('M_{0F}'); ylabel('PDF Value'); hold on; vline(M0F_Prediction, 'r', 'Pred')
subplot(2,4,3); plot(M0S_Vector,exp(P_M0S),'k'); xlabel('M_{0S}'); ylabel('PDF Value'); hold on; vline(M0S_Prediction, 'r', 'Pred')
subplot(2,4,5); plot(kFS_Vector,exp(P_kFS),'k'); xlabel('k_{FS}'); ylabel('PDF Value'); hold on; vline(kFS_Prediction, 'r', 'Pred')
subplot(2,4,6); plot(kSF_Vector,exp(P_kSF),'k'); xlabel('k_{SF}'); ylabel('PDF Value'); hold on; vline(kSF_Prediction, 'r', 'Pred')
subplot(2,4,7); plot(T2S_Vector,exp(P_T2S),'k'); xlabel('T_{2S}'); ylabel('PDF Value'); hold on; vline(T2_S, 'r', 'Pred')
subplot(2,4,8); plot(T2F_Vector,exp(P_T2F),'k'); xlabel('T_{2F}'); ylabel('PDF Value'); hold on; vline(T2_F, 'r', 'Pred')

%% Direct sampling of cost-function in 2D.

Steps = 500;

T1S_Vector = linspace(Lower(1),Upper(1),Steps);
T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps);
kFS_Vector = linspace(Lower(5),Upper(5),Steps);
kSF_Vector = linspace(Lower(6),Upper(6),Steps);
T2S_Vector = linspace(Lower(7),Upper(7),Steps);
T2F_Vector = linspace(Lower(8),Upper(8),Steps);

% T1S_Vector = linspace(0.2,0.3,Steps);
% T1F_Vector = linspace(0.15,0.25,Steps);
% M0F_Vector = linspace(0.1,0.2,Steps);
% %M0S_Vector = linspace(0.1,0.15,Steps);
% kFS_Vector = linspace(10,11,Steps);
% kSF_Vector = linspace(4.5,5,Steps);
% T2S_Vector = linspace(0.05,0.1,Steps);
% T2F_Vector = linspace(0.01,0.05,Steps);

P_T1S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T1F = zeros(length(T1S_Vector),length(T1S_Vector));
%P_M0S = zeros(length(T1S_Vector),length(T1S_Vector));
P_kFS = zeros(length(T1S_Vector),length(T1S_Vector));
P_kSF = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2S = zeros(length(T1S_Vector),length(T1S_Vector));
P_T2F = zeros(length(T1S_Vector),length(T1S_Vector));
%P_T1 = zeros(length(T1S_Vector),length(T1S_Vector));
%P_k = zeros(length(T1S_Vector),length(T1S_Vector));
%P_T2 = zeros(length(T1S_Vector),length(T1S_Vector));

disp (['Estimated Run Time: ' num2str(0.0105*Steps^2) ' seconds']);
tic
for ii = 1:length(M0F_Vector)
    tic
    for jj = 1:length(T1S_Vector)
        
        P_T1S(ii,jj) = (logpdf([T1S_Vector(ii) T1F_Prediction M0F_Vector(jj) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T1F(ii,jj) = (logpdf([T1S_Prediction T1F_Vector(ii) M0F_Vector(jj) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_M0S(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(jj) M0S_Vector(ii) kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kFS(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Vector(jj) kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_kSF(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2S(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2S_Vector(jj) T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        P_T2F(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Vector(ii) M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T1(ii,jj) = (logpdf([T1S_Vector(ii) T1F_Vector(jj) M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_k(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Vector(ii) kSF_Vector(jj) T2_S T2_F], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        %P_T2(ii,jj) = (logpdf([T1S_Prediction T1F_Prediction M0F_Prediction M0S_Prediction kFS_Prediction kSF_Prediction T2S_Vector(ii) T2F_Vector(jj)], Lower, Upper, FA_SPGR, FA_SSFP, TR_SPGR, TR_SSFP, Data, Sigma));
        
    end
   toc 
end
toc

Matrix_kFS = exp(P_kFS); Matrix_kSF = exp(P_kSF);
Matrix_T1F = exp(P_T1F); Matrix_T1S = exp(P_T1S);
Matrix_T2F = exp(P_T2F); Matrix_T2S = exp(P_T2S);
Threshold = 0.99999;

[RowX_kFS, ColX_kFS] = find(Matrix_kFS > Threshold);
[RowX_kSF, ColX_kSF] = find(Matrix_kSF > Threshold);
[RowX_T1F, ColX_T1F] = find(Matrix_T1F > Threshold);
[RowX_T1S, ColX_T1S] = find(Matrix_T1S > Threshold);
[RowX_T2F, ColX_T2F] = find(Matrix_T2F > Threshold);
[RowX_T2S, ColX_T2S] = find(Matrix_T2S > Threshold);

figure(5)
subplot(3,2,1)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1S)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1S}'); hold on; plot(M0_F,T1_S,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1S_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1S), T1S_Vector(RowX_T1S), 'rx')

subplot(3,2,2)
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],exp(P_T1F)); colorbar; shading interp; xlabel('M_{0F}'); ylabel('T_{1F}'); hold on; plot(M0_F,T1_F,'w.', 'MarkerSize', 20); plot(M0F_Prediction,T1F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(M0F_Vector(ColX_T1F), T1F_Vector(RowX_T1F), 'rx')

% subplot(5,2,3)
% imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_M0S)); colorbar; shading interp; xlabel('M_{0S}'); ylabel('M_{0F}')

subplot(3,2,3)
imagesc([min(kFS_Vector) max(kFS_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kFS)); colorbar; shading interp; xlabel('k_{FS}'); ylabel('M_{0F}'); hold on; plot(k_FS,M0_F,'w.', 'MarkerSize', 20); plot(kFS_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(kFS_Vector(ColX_kFS), M0F_Vector(RowX_kFS), 'rx')

subplot(3,2,4)
imagesc([min(kSF_Vector) max(kSF_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_kSF)); colorbar; shading interp; xlabel('k_{SF}'); ylabel('M_{0F}'); hold on; plot(k_SF,M0_F,'w.', 'MarkerSize', 20); plot(kSF_Prediction,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(kSF_Vector(ColX_kSF), M0F_Vector(RowX_kSF), 'rx')

subplot(3,2,5)
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2S)); colorbar; shading interp; xlabel('T_{2S}'); ylabel('M_{0F}'); hold on; plot(T2_S,M0_F,'w.', 'MarkerSize', 20); plot(T2_S,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(T2S_Vector(ColX_T2S), M0F_Vector(RowX_T2S), 'rx')

subplot(3,2,6)
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],exp(P_T2F)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('M_{0F}'); hold on; plot(T2_F,M0_F,'w.', 'MarkerSize', 20); plot(T2_F,M0F_Prediction,'k.', 'MarkerSize', 20);
hold on; scatter(T2F_Vector(ColX_T2F), M0F_Vector(RowX_T2F), 'rx')

% subplot(5,2,8)
% imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],exp(P_T1)); colorbar; shading interp; xlabel('T_{1F}'); ylabel('T_{1S}')
% subplot(5,2,9)
% imagesc([min(kSF_Vector) max(kSF_Vector)],[min(kFS_Vector) max(kFS_Vector)],exp(P_k)); colorbar; shading interp; xlabel('k_{SF}'); ylabel('k_{FS}')
% subplot(5,2,10)
% imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],exp(P_T2)); colorbar; shading interp; xlabel('T_{2F}'); ylabel('T_{2S}')

%% Donald's suggestion.
% 
% figure(5); hold on
% % subplot(5,2,1)
% % S1 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S1 = 0.01*S1+repmat(Start(1:2),size(S1,1),1); S1 = [ S1 repmat(Start(3:end),size(S1,1),1) ];
% % for ii=1:size(S1,1),lnP_val(ii) = lnP(S1(ii,:)); end
% % scatter(S1(:,1), S1(:,2), 10, exp(lnP_val), 'filled'); shading interp; colorbar; %T1s.
% % xlim([min(S1(:,1)) max(S1(:,1))]); ylim([min(S1(:,2)) max(S1(:,2))]); xlabel('T_{1S}'); ylabel('T_{1F}');
% % 
% % subplot(5,2,2)
% % S2 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S2 = 1*S2+repmat(Start(5:6),size(S2,1),1); S2 = [ repmat(Start(1:4),size(S2,1),1) S2  repmat(Start(7:end),size(S2,1),1)];
% % for ii=1:size(S2,1),lnP_val(ii) = lnP(S2(ii,:)); end
% % scatter (S2(:,5), S2(:,6), 50, exp(lnP_val), '.'); shading interp; colorbar; title('ks') %ks.
% % xlim([min(S2(:,5)) max(S2(:,5))]); ylim([min(S2(:,6)) max(S1(:,6))])
% % 
% % subplot(5,2,3)
% % S3 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S3 = 0.01*S3+repmat(Start(7:8),size(S3,1),1); S3 = [ repmat(Start(1:6),size(S3,1),1) S3];
% % for ii=1:size(S3,1),lnP_val(ii) = lnP(S3(ii,:)); end
% % scatter (S3(:,7), S3(:,8), 50, exp(lnP_val), '.'); shading interp; colorbar; title('T_{2}s') %T2s.
% % xlim([min(S3(:,7)) max(S1(:,7))]); ylim([min(S3(:,8)) max(S3(:,8))])
% 
% % Examination of M0F correlations to all other parameters.
% 
% subplot(3,2,1)
% S4 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S4 = 0.01*S4+repmat(Start([1 3]),size(S4,1),1); S4 = [ S4(:,1) repmat(Start(2),size(S4,1),1) S4(:,2) repmat(Start(4:end),size(S4,1),1) ];
% for ii=1:size(S4,1),lnP_val(ii) = lnP(S4(ii,:)); end
% scatter (S4(:,1), S4(:,3), 50, exp(lnP_val), '.'); shading interp; colorbar; title('T_{1S}-M_{0F}')
% xlim([min(S4(:,1)) max(S4(:,1))]); ylim([min(S4(:,3)) max(S4(:,3))])
% 
% subplot(3,2,2)
% S5 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S5 = 0.01*S5+repmat(Start([2 3]),size(S5,1),1); S5 = [ repmat(Start(1),size(S5,1),1) S5 repmat(Start(4:end),size(S5,1),1) ];
% for ii=1:size(S5,1),lnP_val(ii) = lnP(S5(ii,:)); end
% scatter (S5(:,2), S5(:,3), 50, exp(lnP_val), '.'); shading interp; colorbar; title('T_{1F}-M_{0F}')
% xlim([min(S5(:,2)) max(S5(:,2))]); ylim([min(S5(:,3)) max(S5(:,3))])
% 
% % subplot(5,2,6)
% % S6 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S6 = 0.01*S6+repmat(Start([3 4]),size(S6,1),1); S6 = [ repmat(Start(1:2),size(S6,1),1) S6 repmat(Start(5:end),size(S6,1),1) ];
% % for ii=1:size(S6,1),lnP_val(ii) = lnP(S6(ii,:)); end
% % scatter (S6(:,3), S6(:,4), 50, exp(lnP_val), '.'); shading interp; colorbar; title('M_{0F}-M_{0S}')
% % xlim([min(S6(:,3)) max(S6(:,3))]); ylim([min(S6(:,4)) max(S6(:,4))])
% 
% subplot(3,2,3)
% S7 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S7 = 0.01*S7+repmat(Start([3 5]),size(S7,1),1); S7 = [ repmat(Start(1:2),size(S7,1),1) S7(:,1) repmat(Start(4),size(S7,1),1) S7(:,2) repmat(Start(6:end),size(S7,1),1) ];
% for ii=1:size(S7,1),lnP_val(ii) = lnP(S7(ii,:)); end
% scatter (S7(:,3), S7(:,5), 50, exp(lnP_val), '.'); shading interp; colorbar; title('M_{0F}-k_{FS}')
% xlim([min(S7(:,3)) max(S7(:,3))]); ylim([min(S7(:,5)) max(S7(:,5))])
% 
% subplot(3,2,4)
% S8 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S8 = 0.01*S8+repmat(Start([3 6]),size(S8,1),1); S8 = [ repmat(Start(1:2),size(S8,1),1) S8(:,1) repmat(Start(4:5),size(S8,1),1) S8(:,2) repmat(Start(7:end),size(S8,1),1) ];
% for ii=1:size(S8,1),lnP_val(ii) = lnP(S8(ii,:)); end
% scatter (S8(:,3), S8(:,6), 50, exp(lnP_val), '.'); shading interp; colorbar; title('M_{0F}-k_{SF}')
% xlim([min(S8(:,3)) max(S8(:,3))]); ylim([min(S8(:,6)) max(S8(:,6))])
% 
% subplot(3,2,5)
% S9 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S9 = 0.01*S9+repmat(Start([3 7]),size(S9,1),1); S9 = [ repmat(Start(1:2),size(S9,1),1) S9(:,1) repmat(Start(4:6),size(S9,1),1) S9(:,2) repmat(Start(8),size(S9,1),1) ];
% for ii=1:size(S9,1),lnP_val(ii) = lnP(S9(ii,:)); end
% scatter (S9(:,3), S9(:,7), 50, exp(lnP_val), '.'); shading interp; colorbar; title('M_{0F}-T_{2S}')
% xlim([min(S9(:,3)) max(S9(:,3))]); ylim([min(S9(:,7)) max(S9(:,7))])
% 
% subplot(3,2,6)
% S10 = [ reshape(repmat(-1:0.02:1, 101,1),[],1) reshape(repmat((-1:0.02:1)',1,101),[], 1) ]; S10 = 0.01*S10+repmat(Start([3 8]),size(S10,1),1); S10 = [ repmat(Start(1:2),size(S10,1),1) S10(:,1) repmat(Start(4:7),size(S10,1),1) S10(:,2) ];
% for ii=1:size(S10,1),lnP_val(ii) = lnP(S10(ii,:)); end
% scatter (S10(:,3), S10(:,8), 50, exp(lnP_val), '.'); shading interp; colorbar; title('M_{0F}-T_{2F}')
% xlim([min(S10(:,3)) max(S10(:,3))]); ylim([min(S10(:,8)) max(S10(:,8))])

