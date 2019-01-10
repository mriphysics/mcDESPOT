%% Generates heat maps to depict cost-function space. (MaxAngle is in degrees and phase-cycling is in radians).

close all; clear all;

% Tissue and sequence parameters.
%TR_SPGR = 5e-3; TR_SSFP = 5e-3; FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([10 13 17 20 23 30 43 60]); FA_SSFP180 = deg2rad([10 13 17 20 23 30 43 60]); MaxAngle = 60;
TR_SPGR = 7e-3; TR_SSFP = 7e-3; FA_SPGR = deg2rad([6 8 10 12 14 16]); FA_SSFP180 = deg2rad([15 25 35 45 55 65]); FA_SSFP0 = deg2rad([25 55]); MaxAngle = 65;

T1_S = 0.965; T1_F = 0.465; T1_B = 1; T2_S = 0.09; T2_F = 0.012; M0_F = 0.2; k_FS = 8; M0_S = 0.7; M0_B = 0.1; k_SB = 5; k_FB = 5; k_SF = (M0_F*k_FS)/M0_S;
Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

% TR_SPGR = 5.2e-3; TR_SSFP = 5.2e-3; MaxAngle = 50;
% FA_SPGR = deg2rad([3 4 5 6 7 9 13 18]); FA_SSFP0 = deg2rad([2 5 10 15 20 30 40 50]); FA_SSFP180 = deg2rad([2 5 10 15 20 30 40 50]);
% T1_S = 1.15; T1_F = 0.4; T1_B = 1; T2_S = 0.08; T2_F = 0.02; M0_F = 0.25; M0_S = 0.55; M0_B = 0.2; k_FS = 9; k_SB = 5; k_FB = 5;
% PC1 = 0; PC2 = pi; Sigma = 1/300; %SNR = 30.

%% Specify ground-truth signals.

% Normalise, concatenate and change between CS and US (B1 or TRF).
SPGR_Data = B1_SPGR_SteadyState(FA_SPGR,TR_SPGR,'T1_S',T1_S,'T1_F',T1_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_0 = B1_SSFP_SteadyState(FA_SSFP0,TR_SSFP,PC1,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data_180 = B1_SSFP_SteadyState(FA_SSFP180,TR_SSFP,PC2,'T1_S',T1_S,'T2_S',T2_S,'T1_F',T1_F,'T2_F',T2_F,'T1_B',T1_B,'M0_B',M0_B,'M0_F',M0_F,'M0_S',M0_S,'k_FS',k_FS,'k_SB',k_SB,'k_FB',k_FB);
SSFP_Data = [SSFP_Data_0 ; SSFP_Data_180];
Data = [SPGR_Data; SSFP_Data];

SNR = 30; Sigma = mean(SPGR_Data)/SNR;

%% Direct sampling of cost-function in 2D.

Upper = [1.25 0.5 0.30 0.80 20 10 0.15 0.03 pi]; Lower = [0.20 0.1 0.05 0.05 1 1 0.04 0.01 -pi]; % Narrower.
%Upper = [3 1.5 1 1 40 20 0.2 0.05 pi]; Lower = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 -pi]; % Wider.
Steps = 500; nTrials = 1;

T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
M0F_Vector = linspace(Lower(3),Upper(3),Steps); M0S_Vector = linspace(Lower(4),Upper(4),Steps);
kFS_Vector = linspace(Lower(5),Upper(5),Steps); kSF_Vector = linspace(Lower(6),Upper(6),Steps);
T2S_Vector = linspace(Lower(7),Upper(7),Steps); T2F_Vector = linspace(Lower(8),Upper(8),Steps);
Delta_Vector = linspace(Lower(9),Upper(9),Steps);

P_T1S = zeros(length(T1S_Vector),length(M0F_Vector),nTrials);
P_T1F = zeros(length(T1F_Vector),length(M0F_Vector),nTrials);
P_T2S = zeros(length(M0F_Vector),length(T1S_Vector),nTrials);
P_T2F = zeros(length(M0F_Vector),length(T1S_Vector),nTrials);
P_T1 = zeros(length(T1S_Vector),length(M0F_Vector),nTrials);
P_T2 = zeros(length(T2S_Vector),length(M0F_Vector),nTrials);

for ii = 1:length(T1S_Vector)
    tic
    disp(['Step Number: ', num2str(ii), '/', num2str(Steps), '.'])
    for jj = 1:length(T1F_Vector)
        
        for nn = 1:nTrials
            % Choose either Prediction or GT as first combination.
            if nn == 1
                
                %T1F_Rand = T1F_Prediction; T1S_Rand = T1S_Prediction;
                %T2F_Rand = T2_F; T2S_Rand = T2_S;
                %kFS_Rand = kFS_Prediction; kSF_Rand = kSF_Prediction;
                %M0F_Rand = M0F_Prediction; M0S_Rand = M0S_Prediction;
                T1F_Rand = T1_F; T1S_Rand = T1_S;
                T2F_Rand = T2_F; T2S_Rand = T2_S;
                kFS_Rand = k_FS; kSF_Rand = k_SF;
                M0F_Rand = M0_F; M0S_Rand = M0_S;
                Delta_Rand = Delta;
                 
            else
                
                % Draw random parameter values from a uniform distribution.
                T1S_Rand = (Upper(1) - Lower(1)) .* rand(1,1) + Lower(1);
                T1F_Rand = (Upper(2) - Lower(2)) .* rand(1,1) + Lower(2);
                M0F_Rand = (Upper(3) - Lower(3)) .* rand(1,1) + Lower(3);
                M0S_Rand = (Upper(4) - Lower(4)) .* rand(1,1) + Lower(4);
                kFS_Rand = (Upper(5) - Lower(5)) .* rand(1,1) + Lower(5);
                kSF_Rand = (Upper(6) - Lower(6)) .* rand(1,1) + Lower(6);
                T2S_Rand = (Upper(7) - Lower(7)) .* rand(1,1) + Lower(7);
                T2F_Rand = (Upper(8) - Lower(8)) .* rand(1,1) + Lower(8);
                Delta_Rand = (Upper(9) - Lower(9)) .* rand(1,1) + Lower(9);
                
            end
            
            P_T1S(ii,jj,nn) = (logpdf([T1S_Vector(ii) T1F_Rand M0F_Vector(jj) M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T1F(ii,jj,nn) = (logpdf([T1S_Rand T1F_Vector(ii) M0F_Vector(jj) M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T2S(ii,jj,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Vector(ii) M0S_Rand kFS_Rand kSF_Rand T2S_Vector(jj) T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T2F(ii,jj,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Vector(ii) M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Vector(jj) Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T1(ii,jj,nn) = (logpdf([T1S_Vector(ii) T1F_Vector(jj) M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Rand T2F_Rand Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            P_T2(ii,jj,nn) = (logpdf([T1S_Rand T1F_Rand M0F_Rand M0S_Rand kFS_Rand kSF_Rand T2S_Vector(ii) T2F_Vector(jj) Delta_Rand], Lower, Upper, FA_SPGR, FA_SSFP0, FA_SSFP180, TR_SPGR, TR_SSFP, Data, Sigma));
            
        end
        
    end
    toc
end

%% Heat map analysis.

figure(1); subplot(2,3,1);
MaxCF_T1 = max(exp(P_T1),[],3);
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],MaxCF_T1); colormap(viridis); colorbar; shading interp; xlabel('T_{1F}','FontSize',12); ylabel('T_{1S}','FontSize',12); hold on;
plot(T1_F,T1_S,'w.', 'MarkerSize', 20); %plot(T1F_Prediction,T1S_Prediction,'m.','MarkerSize',20);
% [Row, Col] = find(ismember(MaxCF_T1, max(MaxCF_T1(:))));
% Max_T1F = T1F_Vector(Col); Max_T1S = T1S_Vector(Row); plot(Max_T1F,Max_T1S,'ro','LineWidth',2);

subplot(2,3,2);
MaxCF_T2 = max(exp(P_T2),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],MaxCF_T2); colormap(viridis); colorbar; shading interp; xlabel('T_{2F}','FontSize',12); ylabel('T_{2S}','FontSize',12); hold on;
plot(T2_F,T2_S,'w.', 'MarkerSize', 20);
% [Row, Col] = find(ismember(MaxCF_T2, max(MaxCF_T2(:))));
% Max_T2F = T2F_Vector(Col); Max_T2S = T2S_Vector(Row); plot(Max_T2F,Max_T2S,'ro','LineWidth',2);

subplot(2,3,3);
MaxCF_T1F = max(exp(P_T1F),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],MaxCF_T1F); colormap(viridis); colorbar; shading interp; xlabel('M_{0F}','FontSize',12); ylabel('T_{1F}','FontSize',12); hold on;
plot(M0_F,T1_F,'w.', 'MarkerSize', 20); %plot(M0F_Prediction,T1F_Prediction,'m.','MarkerSize',20);
% [Row, Col] = find(ismember(MaxCF_T1F, max(MaxCF_T1F(:))));
% Max_M0F = M0F_Vector(Col); Max_T1F = T1F_Vector(Row); plot(Max_M0F,Max_T1F,'ro','LineWidth',2);

subplot(2,3,4);
MaxCF_T1S = max(exp(P_T1S),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],MaxCF_T1S); colormap(viridis); colorbar; shading interp; xlabel('M_{0F}','FontSize',12); ylabel('T_{1S}','FontSize',12); hold on;
plot(M0_F,T1_S,'w.', 'MarkerSize', 20); %plot(M0F_Prediction,T1S_Prediction,'m.','MarkerSize',20);
% [Row, Col] = find(ismember(MaxCF_T1S, max(MaxCF_T1S(:))));
% Max_M0F = M0F_Vector(Col); Max_T1S = T1S_Vector(Row); plot(Max_M0F,Max_T1S,'ro','LineWidth',2);

subplot(2,3,5);
MaxCF_T2F = max(exp(P_T2F),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_T2F); colormap(viridis); colorbar; shading interp; xlabel('T_{2F}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on;
plot(T2_F,M0_F,'w.', 'MarkerSize', 20); %plot(T2_F,M0F_Prediction,'m.','MarkerSize',20);
% [Row, Col] = find(ismember(MaxCF_T2F, max(MaxCF_T2F(:))));
% Max_T2F = T2F_Vector(Col); Max_M0F = M0F_Vector(Row); plot(Max_T2F,Max_M0F,'ro','LineWidth',2);

subplot(2,3,6);
MaxCF_T2S = max(exp(P_T2S),[],3);
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_T2S); colormap(viridis); colorbar; shading interp; xlabel('T_{2S}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on;
plot(T2_S,M0_F,'w.', 'MarkerSize', 20); %plot(T2_S,M0F_Prediction,'m.','MarkerSize',20);
% [Row, Col] = find(ismember(MaxCF_T2S, max(MaxCF_T2S(:))));
% Max_T2S = T2S_Vector(Col); Max_M0F = M0F_Vector(Row); plot(Max_T2S,Max_M0F,'ro','LineWidth',2);

%% Outtakes.

% figure(6);
% MaxCF_M0 = max(exp(P_M0),[],3); MeanCF_M0 = mean(exp(P_M0),3);
% subplot(2,1,1); 
% imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_M0); colormap(viridis); colorbar; shading interp; xlabel('M_{0S}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on; plot(M0_S,M0_F,'w.', 'MarkerSize', 20); %plot(M0S_Prediction,M0F_Prediction,'m.', 'MarkerSize', 20);
%[M0S_Mesh, M0F_Mesh] = meshgrid(M0S_Vector, M0F_Vector);
%[C,h] = contour(M0S_Vector, M0F_Vector, (M0F_Mesh./(M0F_Mesh + M0S_Mesh)),'LineStyle','--','LineWidth',2,'ShowText','on','LineColor','w');
%clabel(C,h,'FontSize',12,'Color','w')
%MWF_Prediction = M0F_Prediction/(M0F_Prediction+M0S_Prediction); M0S_Test = M0S_Vector; M0F_Test = (1/((1/MWF_Prediction)-1))*M0S_Test;
%line(M0S_Test,M0F_Test,'LineWidth',2,'LineStyle','--','Color','m')
%[Row, Col] = find(ismember(MaxCF_M0, max(MaxCF_M0(:))));
%Max_M0S = M0S_Vector(Col); Max_M0F = M0F_Vector(Row);
%plot(Max_M0S,Max_M0F,'ro','LineWidth',2);

% subplot(2,1,2); imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MaxCF_M0); colormap(viridis); colorbar; shading interp; xlabel('M_{0S}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on; plot(M0_S,M0_F,'w.', 'MarkerSize', 20); %plot(M0S_Prediction,M0F_Prediction,'m.', 'MarkerSize', 20);
% [M0S_Mesh, M0F_Mesh] = meshgrid(M0S_Vector, M0F_Vector);
% [C,h] = contour(M0S_Vector, M0F_Vector, (M0F_Mesh./(M0F_Mesh + M0S_Mesh)),'LineStyle','--','LineWidth',2,'ShowText','on','LineColor','w');
% clabel(C,h,'FontSize',12,'Color','w')
% MWF_Prediction = M0F_Prediction/(M0F_Prediction+M0S_Prediction); M0S_Test = M0S_Vector; M0F_Test = (1/((1/MWF_Prediction)-1))*M0S_Test;
% line(M0S_Test,M0F_Test,'LineWidth',2,'LineStyle','--','Color','m')
% [Row, Col] = find(ismember(MaxCF_M0, max(MaxCF_M0(:))));
% Max_M0S = M0S_Vector(Col); Max_M0F = M0F_Vector(Row);
% plot(Max_M0S,Max_M0F,'ro','LineWidth',2);

%figure(7); 
%MinSSE_M0 = min(SSE_M0,[],3);
% G_Factor = 0.1; MinSSE_M0 = imgaussfilt(MaxCF_M0,G_Factor);
%imagesc([min(M0S_Vector) max(M0S_Vector)],[min(M0F_Vector) max(M0F_Vector)],MinSSE_M0); colormap(magma); colorbar; shading interp; xlabel('M_{0S}','FontSize',12); ylabel('M_{0F}','FontSize',12); hold on; plot(M0_S,M0_F,'w.', 'MarkerSize', 20); %plot(M0S_Prediction,M0F_Prediction,'m.', 'MarkerSize', 20);

%end

%figure(10); G_Factor = 2.5;
% MaxCF_M0_Filtered = imgaussfilt(MaxCF_M0,G_Factor);
% while MaxCF_M0_Filtered(:) < 0.5*MaxVal_M0
%     G_Factor = G_Factor - 0.1;
%     MaxCF_M0_Filtered = imgaussfilt(MaxCF_M0,G_Factor);
% end

% textLabel1 = sprintf('(%.2f, %.2f)', Max_M0S, Max_M0F); text(Max_M0S-0.04, Max_M0F-0.005, 'True Max','Color','r','FontSize',12);
% textLabel2 = sprintf('(%.2f, %.2f)', M0_S, M0_F); text(M0_S-0.025, M0_F-0.005, '3P GT','Color','w','FontSize',12);
% textLabel3 = sprintf('(%.2f, %.2f)', M0S_Prediction, M0F_Prediction); text(M0S_Prediction-0.04, M0F_Prediction-0.005, 'Prediction','Color','w','FontSize',12);
% MWF_Test = M0F_Prediction/(M0F_Prediction+M0S_Prediction); M0F_Test = (1/((1/MWF_Test)-1))*M0S_Vector;
% hold on; line(M0S_Vector,M0F_Test,'Color','red','LineStyle','--','LineWidth',2)

% figure(1); subplot(2,2,3);
% MeanCF_EX = mean(exp(P_EX),3); SDCF_EX = std(exp(P_EX),0,3); MaxCF_EX = max(exp(P_EX),[],3);
% MaxVal_EX = max(max(MaxCF_EX)); [Values_EX, Index_EX] = find(MaxCF_EX > 0.5); Percentage_EX = (length(Values_EX)/(size(P_EX,1)*size(P_EX,2)))*100;
% imagesc([min(kSF_Vector) max(kSF_Vector)],[min(kFS_Vector) max(kFS_Vector)],MaxCF_EX); colormap(viridis); colorbar; shading interp; xlabel('k_{SF}','FontSize',12); ylabel('k_{FS}','FontSize',12); hold on; plot(k_SF,k_FS,'w.', 'MarkerSize', 20); plot(kSF_Prediction,kFS_Prediction,'m.', 'MarkerSize', 20);
