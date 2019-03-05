%%% Generates search space heat maps in Figures 3-4. %%%

acquisition = 'Bouhrara';
exchange = 'off';

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

Steps_MP = 49; nTrials_MP = 100000; Steps = 499; nTrials = 1;

%% Maximum projection heat maps.

PointC = [0 0.7 0]; col_axis = [0 1];

switch exchange
    case 'on'
        Upper = [2 0.7 0.3 16 0.15 0.03]; Lower = [0.8 0.2 0 0 0.03 0];
        T1S_Vector = linspace(Lower(1),Upper(1),Steps_MP); T1F_Vector = linspace(Lower(2),Upper(2),Steps_MP);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps_MP); 
        T2S_Vector = linspace(Lower(5),Upper(5),Steps_MP); T2F_Vector = linspace(Lower(6),Upper(6),Steps_MP);
    case 'off'
        Upper = [2 0.7 0.3 0.15 0.03]; Lower = [0.8 0.2 0 0.03 0];
        T1S_Vector = linspace(Lower(1),Upper(1),Steps_MP); T1F_Vector = linspace(Lower(2),Upper(2),Steps_MP);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps_MP);
        T2S_Vector = linspace(Lower(4),Upper(4),Steps_MP); T2F_Vector = linspace(Lower(5),Upper(5),Steps_MP);
end

figure(1); subplot(3,6,1); %MaxCF_T1_MP = max(exp(P_T1),[],3);
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{1F} (s)','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T1_F,T1_S,'.','Color',PointC,'MarkerSize', 20);

subplot(3,6,2); %MaxCF_T2_MP = max(exp(P_T2),[],3);
MaxCF_T2_MP_NoiseInd(isnan(MaxCF_T2_MP_NoiseInd)) = Inf;
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],(MaxCF_T2_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',16); ylabel('T_{2S} (s)','FontSize',16); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,T2_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,3); %MaxCF_T1F_MP = max(exp(P_T1F),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],(MaxCF_T1F_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1F} (s)','FontSize',16); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,4); %MaxCF_T1S_MP = max(exp(P_T1S),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1S_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,5); %MaxCF_T2F_MP = max(exp(P_T2F),[],3);
MaxCF_T2F_MP_NoiseInd(isnan(MaxCF_T2F_MP_NoiseInd)) = Inf;
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2F_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',14); ylabel('MWF','FontSize',14); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,M0_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,6); %MaxCF_T2S_MP = max(exp(P_T2S),[],3);
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2S_MP_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2S} (s)','FontSize',16); ylabel('MWF','FontSize',16); hold on;
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_S,M0_F,'.','Color',PointC, 'MarkerSize', 20);

%% GT heat maps.

PointC = [0 0.7 0];

switch exchange
    case 'on'
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps);
        T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);
    case 'off'
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps);
        T2S_Vector = linspace(Lower(4),Upper(4),Steps); T2F_Vector = linspace(Lower(5),Upper(5),Steps);
end

T1S_LB = Bounds_SNR100(:,1); T1S_UB = Bounds_SNR100(:,2); T1F_LB = Bounds_SNR100(:,3); T1F_UB = Bounds_SNR100(:,4);
T2S_LB = Bounds_SNR100(:,5); T2S_UB = Bounds_SNR100(:,6); T2F_LB = Bounds_SNR100(:,7); T2F_UB = Bounds_SNR100(:,8); 
M0F_LB = Bounds_SNR100(:,9); M0F_UB = Bounds_SNR100(:,10);

figure(1); subplot(3,6,7); %MaxCF_T1_GT = max(exp(P_T1_GT),[],3);
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{1F} (s)','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T1F_LB(rr) T1S_LB(rr) (T1F_UB(rr)-T1F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T1_F,T1_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,8); %MaxCF_T2_GT = max(exp(P_T2_GT),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],(MaxCF_T2_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',16); ylabel('T_{2S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2F_LB(rr) T2S_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (T2S_UB(rr)-T2S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,T2_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,9); %MaxCF_T1F_GT = max(exp(P_T1F_GT),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],(MaxCF_T1F_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1F} (s)','FontSize',14); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[M0F_LB(rr) T1F_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1F_UB(rr)-T1F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,10); %MaxCF_T1S_GT = max(exp(P_T1S_GT),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1S_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[M0F_LB(rr) T1S_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,11); %MaxCF_T2F_GT = max(exp(P_T2F_GT),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2F_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',16); ylabel('MWF','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2F_LB(rr) M0F_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,M0_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,12); %MaxCF_T2S_GT = max(exp(P_T2S_GT),[],3);
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2S_GT_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2S} (s)','FontSize',16); ylabel('MWF','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2S_LB(rr) M0F_LB(rr) (T2S_UB(rr)-T2S_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_S,M0_F,'.','Color',PointC, 'MarkerSize', 20);

%% SRC heat maps.

PointC = [0 0.7 0];

switch exchange
    case 'on'
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps);
        T2S_Vector = linspace(Lower(5),Upper(5),Steps); T2F_Vector = linspace(Lower(6),Upper(6),Steps);
    case 'off'
        T1S_Vector = linspace(Lower(1),Upper(1),Steps); T1F_Vector = linspace(Lower(2),Upper(2),Steps);
        M0F_Vector = linspace(Lower(3),Upper(3),Steps);
        T2S_Vector = linspace(Lower(4),Upper(4),Steps); T2F_Vector = linspace(Lower(5),Upper(5),Steps);
end

T1S_LB = Bounds_SNR100(:,1); T1S_UB = Bounds_SNR100(:,2); T1F_LB = Bounds_SNR100(:,3); T1F_UB = Bounds_SNR100(:,4);
T2S_LB = Bounds_SNR100(:,5); T2S_UB = Bounds_SNR100(:,6); T2F_LB = Bounds_SNR100(:,7); T2F_UB = Bounds_SNR100(:,8); 
M0F_LB = Bounds_SNR100(:,9); M0F_UB = Bounds_SNR100(:,10);

figure(1); subplot(3,6,13); %MaxCF_T1_SRC = max(exp(P_T1_SRC),[],3);
imagesc([min(T1F_Vector) max(T1F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1_SRC_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{1F} (s)','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T1F_LB(rr) T1S_LB(rr) (T1F_UB(rr)-T1F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T1_F,T1_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,14); %MaxCF_T2_SRC = max(exp(P_T2_SRC),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(T2S_Vector) max(T2S_Vector)],(MaxCF_T2_SRC_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',16); ylabel('T_{2S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2F_LB(rr) T2S_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (T2S_UB(rr)-T2S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,T2_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,15); %MaxCF_T1F_SRC = max(exp(P_T1F_SRC),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1F_Vector) max(T1F_Vector)],(MaxCF_T1F_SRC_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1F} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[M0F_LB(rr) T1F_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1F_UB(rr)-T1F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,16); %MaxCF_T1S_SRC = max(exp(P_T1S_SRC),[],3);
imagesc([min(M0F_Vector) max(M0F_Vector)],[min(T1S_Vector) max(T1S_Vector)],(MaxCF_T1S_SRC_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('MWF','FontSize',16); ylabel('T_{1S} (s)','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[M0F_LB(rr) T1S_LB(rr) (M0F_UB(rr)-M0F_LB(rr)) (T1S_UB(rr)-T1S_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(M0_F,T1_S,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,17); %MaxCF_T2F_SRC = max(exp(P_T2F_SRC),[],3);
imagesc([min(T2F_Vector) max(T2F_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2F_SRC_NoiseInd),col_axis); colormap(flipud(magma)); shading interp; xlabel('T_{2F} (s)','FontSize',16); ylabel('MWF','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2F_LB(rr) M0F_LB(rr) (T2F_UB(rr)-T2F_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_F,M0_F,'.','Color',PointC, 'MarkerSize', 20);

subplot(3,6,18); %MaxCF_T2S_SRC = max(exp(P_T2S_SRC),[],3);
imagesc([min(T2S_Vector) max(T2S_Vector)],[min(M0F_Vector) max(M0F_Vector)],(MaxCF_T2S_SRC_NoiseInd),col_axis); colormap(flipud(magma)); colorbar; shading interp; xlabel('T_{2S} (s)','FontSize',16); ylabel('MWF','FontSize',16); hold on;
for rr = 1:length(T1S_LB)
    rectangle('Position',[T2S_LB(rr) M0F_LB(rr) (T2S_UB(rr)-T2S_LB(rr)) (M0F_UB(rr)-M0F_LB(rr))],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1)
end
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
plot(T2_S,M0_F,'.','Color',PointC, 'MarkerSize', 20);

colormap(flipud(magma)); 
hcb = colorbar('Position',[0.915318230852205 0.111888111888112 0.0124056094929895 0.813519813519814],'FontSize',16); hcb.FontSize = 18;
%hcb = colorbar('Position',[0.915318230852205,0.41609977324263,0.01240560949299,0.509308152165296],'FontSize',16); hcb.FontSize = 18;
colorTitleHandle = get(hcb,'Title'); titleString = '\eta';
set(colorTitleHandle,'String',titleString,'FontSize',30); 
