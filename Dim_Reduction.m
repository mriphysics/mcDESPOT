%%% Analysis of mcDESPOT search space using dimensionality reduction. 
%%% Requires toolbox from Laurens van der Maaten (https://lvdmaaten.github.io/drtoolbox/).

% Ground-truth.
T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; % Two-pool.
T1 = 1; T2 = 0.1; M0 = 1; % Single-pool.

% Modify data to include ground-truth. Matrices are generated using
% mcDESPOT_SP.m for two-pool or mcDESPOT_SP_SinglePool.m for single-pool.
E_kPCA_Matrix_Mod = [E_kPCA_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, k_FS]]; E_Values_Top_Mod = [E_Values_Top ; 1];
NE_kPCA_Matrix_Mod = [NE_kPCA_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, 0]]; NE_Values_Top_Mod = [NE_Values_Top ; 1];
SP_kPCA_Matrix_Mod = [SP_kPCA_Matrix ; [T1, T2, M0]]; SP_Values_Top_Mod = [SP_Values_Top ; 1];

% Perform dimensionality reduction.
[E_MappedData, E_Mapping] = compute_mapping(E_kPCA_Matrix_Mod(:,1:5), 'KernelPCA', 3);
[NE_MappedData, NE_Mapping] = compute_mapping(NE_kPCA_Matrix_Mod(:,1:5), 'KernelPCA', 3);

%% Plot results for model comparison.

az = 1.590000000000001e+02; el = 34.800000000000033;

figure(1)
subplot(1,3,1)
colormap(magma); h = plot3k(SP_kPCA_Matrix_Mod,'ColorData',SP_Values_Top_Mod,'ColorRange',[0.5 1],'Marker',{'o',5},'ColorBar',false,'Labels',{'','T_{1} (s)','T_{2} (s)','M_{0}',''}); tt = title('Single-Pool'); tt.FontSize = 16; axis square; grid on;
get (h);
hold on; 
plot3(SP_kPCA_Matrix_Mod(:,1),SP_kPCA_Matrix_Mod(:,2),ones(size(SP_kPCA_Matrix_Mod(:,3))).* min(SP_kPCA_Matrix_Mod(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(SP_kPCA_Matrix_Mod(:,1))).* min(SP_kPCA_Matrix_Mod(:,1)),SP_kPCA_Matrix_Mod(:,2),SP_kPCA_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(SP_kPCA_Matrix_Mod(:,1),ones(size(SP_kPCA_Matrix_Mod(:,2))).* min(SP_kPCA_Matrix_Mod(:,2)),SP_kPCA_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(SP_kPCA_Matrix_Mod(1001,1),SP_kPCA_Matrix_Mod(1001,2),SP_kPCA_Matrix_Mod(1001,3),75,[0 0.5 0],'diamond','LineWidth',5)
scatter3(SP_kPCA_Matrix_Mod(1001,1),SP_kPCA_Matrix_Mod(1001,2),min(SP_kPCA_Matrix_Mod(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(SP_kPCA_Matrix_Mod(:,1)),SP_kPCA_Matrix_Mod(1001,2),SP_kPCA_Matrix_Mod(1001,3),75,'k','diamond','LineWidth',5);
scatter3(SP_kPCA_Matrix_Mod(1001,1),min(SP_kPCA_Matrix_Mod(:,2)),SP_kPCA_Matrix_Mod(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
view(az,el);
text(0.02,0.98,'(a)','Units','Normalized','VerticalAlignment','Top','FontSize',16)

subplot(1,3,2)
colormap(magma); plot3k(NE_MappedData,'ColorData',NE_Values_Top_Mod,'ColorRange',[0.5 1],'Marker',{'o',5},'ColorBar',false,'Labels',{'','Dim 1','Dim 2','Dim 3',''}); tt = title('2-Pool (Zero Exchange)'); tt.FontSize = 16; axis square; grid on;
hold on; 
view(az,el)
plot3(NE_MappedData(:,1),NE_MappedData(:,2),ones(size(NE_MappedData(:,3))).* min(NE_MappedData(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(NE_MappedData(:,1))).* min(NE_MappedData(:,1)),NE_MappedData(:,2),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(NE_MappedData(:,1),ones(size(NE_MappedData(:,2))).* min(NE_MappedData(:,2)),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),NE_MappedData(1001,3),75,[0 0.5 0],'diamond','LineWidth',5)
scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),min(NE_MappedData(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(NE_MappedData(:,1)),NE_MappedData(1001,2),NE_MappedData(1001,3),75,'k','diamond','LineWidth',5);
scatter3(NE_MappedData(1001,1),min(NE_MappedData(:,2)),NE_MappedData(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
view(az,el);
text(0.02,0.98,'(b)','Units','Normalized','VerticalAlignment','Top','FontSize',16)

subplot(1,3,3)
colormap(magma); plot3k(E_MappedData,'ColorData',E_Values_Top_Mod,'ColorRange',[0.5 1],'Marker',{'o',5},'ColorBar',false,'Labels',{'','Dim 1','Dim 2','Dim 3',''}); tt = title('2-Pool (Non-Zero Exchange)'); tt.FontSize = 16; axis square; grid on;
hold on; 
plot3(E_MappedData(:,1),E_MappedData(:,2),ones(size(E_MappedData(:,3))).* min(E_MappedData(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10); hold on
plot3(ones(size(E_MappedData(:,1))).* min(E_MappedData(:,1)),E_MappedData(:,2),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(E_MappedData(:,1),ones(size(E_MappedData(:,2))).* min(E_MappedData(:,2)),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(E_MappedData(1001,1),E_MappedData(1001,2),E_MappedData(1001,3),75,[0 0.5 0],'diamond','LineWidth',5);
scatter3(E_MappedData(1001,1),E_MappedData(1001,2),min(E_MappedData(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(E_MappedData(:,1)),E_MappedData(1001,2),E_MappedData(1001,3),75,'k','diamond','LineWidth',5);
scatter3(E_MappedData(1001,1),min(E_MappedData(:,2)),E_MappedData(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
view(az,el);
text(0.02,0.98,'(c)','Units','Normalized','VerticalAlignment','Top','FontSize',16)

hcb = colorbar('Position',[0.915846994535519,0.381118881118881,0.012021857923497,0.452214452214452],'FontSize',14); caxis([0.5 1]);
hcb.Ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; hcb.TickLabels = {'0.5','0.6','0.7','0.8','0.9','1.0'};
colorTitleHandle = get(hcb,'Title'); titleString = 'P (a.u.)';
set(colorTitleHandle,'String',titleString,'FontSize',16,'Position',[-28.500000000000227,299.2006993006992,0]);

