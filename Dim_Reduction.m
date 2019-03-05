%%% Analysis of mcDESPOT search space using dimensionality reduction. 
%%% Requires toolbox from Laurens van der Maaten (https://lvdmaaten.github.io/drtoolbox/).

% Ground-truth.
T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; % Two-pool.
T1 = 1; T2 = 0.1; M0 = 1; % Single-pool.

% Modify data to include ground-truth. Matrices are generated using
% mcDESPOT_SP.m for two-pool or mcDESPOT_SP_SinglePool.m for single-pool.
E_kPCA_Matrix_Mod = [E_kPCA_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, k_FS]]; E_Values_Top_Mod = [E_Values_Top_NoiseInd ; 0];
NE_kPCA_Matrix_Mod = [NE_kPCA_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, 0]]; NE_Values_Top_Mod = [NE_Values_Top_NoiseInd ; 0];
SP_kPCA_Matrix_Mod = [SP_kPCA_Matrix ; [T1, T2, M0]]; SP_Values_Top_Mod = [SP_Values_Top_NoiseInd ; 0];

% Perform dimensionality reduction.
[E_MappedData, E_Mapping] = compute_mapping(E_kPCA_Matrix_Mod(:,1:5), 'KernelPCA', 3);
[NE_MappedData, NE_Mapping] = compute_mapping(NE_kPCA_Matrix_Mod(:,1:5), 'KernelPCA', 3);

%% Plot results for model comparison.

az = 1.590000000000001e+02; el = 34.800000000000033;

figure(1)
colormap(flipud(magma))
subplot(1,3,1)
scatter3(SP_kPCA_Matrix_Mod(:,1),SP_kPCA_Matrix_Mod(:,2),SP_kPCA_Matrix_Mod(:,3),repmat(50,numel(SP_kPCA_Matrix_Mod(:,1)),1),(SP_Values_Top_Mod),'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6); tt = title('Single-Pool'); tt.FontSize = 18; axis square; grid on;
xlabel('T_{1} (s)', 'FontSize', 18); ylabel('T_{2} (s)', 'FontSize', 18); zlabel('M_{0}', 'FontSize', 18);
hold on; caxis([0 1])
plot3(SP_kPCA_Matrix_Mod(:,1),SP_kPCA_Matrix_Mod(:,2),ones(size(SP_kPCA_Matrix_Mod(:,3))).* min(SP_kPCA_Matrix_Mod(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(SP_kPCA_Matrix_Mod(:,1))).* min(SP_kPCA_Matrix_Mod(:,1)),SP_kPCA_Matrix_Mod(:,2),SP_kPCA_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(SP_kPCA_Matrix_Mod(:,1),ones(size(SP_kPCA_Matrix_Mod(:,2))).* min(SP_kPCA_Matrix_Mod(:,2)),SP_kPCA_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(SP_kPCA_Matrix_Mod(1001,1),SP_kPCA_Matrix_Mod(1001,2),SP_kPCA_Matrix_Mod(1001,3),100,[0 0.5 0],'diamond','LineWidth',10)
scatter3(SP_kPCA_Matrix_Mod(1001,1),SP_kPCA_Matrix_Mod(1001,2),min(SP_kPCA_Matrix_Mod(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(SP_kPCA_Matrix_Mod(:,1)),SP_kPCA_Matrix_Mod(1001,2),SP_kPCA_Matrix_Mod(1001,3),75,'k','diamond','LineWidth',5);
scatter3(SP_kPCA_Matrix_Mod(1001,1),min(SP_kPCA_Matrix_Mod(:,2)),SP_kPCA_Matrix_Mod(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18);
view(az,el);
text(0.02,0.98,'(a)','Units','Normalized','VerticalAlignment','Top','FontSize',18)

subplot(1,3,2)
scatter3(NE_MappedData(:,1),NE_MappedData(:,2),NE_MappedData(:,3),repmat(50,numel(NE_MappedData(:,1)),1),(NE_Values_Top_Mod),'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6); tt = title('2-Pool (Zero Exchange)'); tt.FontSize = 18; axis square; grid on;
xlabel('Dim 1', 'FontSize', 18); ylabel('Dim 2', 'FontSize', 18); zlabel('Dim 3', 'FontSize', 18);
hold on; caxis([0 1])
plot3(NE_MappedData(:,1),NE_MappedData(:,2),ones(size(NE_MappedData(:,3))).* min(NE_MappedData(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(NE_MappedData(:,1))).* min(NE_MappedData(:,1)),NE_MappedData(:,2),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(NE_MappedData(:,1),ones(size(NE_MappedData(:,2))).* min(NE_MappedData(:,2)),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),NE_MappedData(1001,3),100,[0 0.5 0],'diamond','LineWidth',10)
scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),min(NE_MappedData(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(NE_MappedData(:,1)),NE_MappedData(1001,2),NE_MappedData(1001,3),75,'k','diamond','LineWidth',5);
scatter3(NE_MappedData(1001,1),min(NE_MappedData(:,2)),NE_MappedData(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18);
view(az,el);
text(0.02,0.98,'(b)','Units','Normalized','VerticalAlignment','Top','FontSize',18)

subplot(1,3,3) 
scatter3(E_MappedData(:,1),E_MappedData(:,2),E_MappedData(:,3),repmat(50,numel(E_MappedData(:,1)),1),(E_Values_Top_Mod),'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6); tt = title('2-Pool (Non-Zero Exchange)'); tt.FontSize = 18; axis square; grid on;
xlabel('Dim 1', 'FontSize', 18); ylabel('Dim 2', 'FontSize', 18); zlabel('Dim 3', 'FontSize', 18);
hold on; caxis([0 1])
plot3(E_MappedData(:,1),E_MappedData(:,2),ones(size(E_MappedData(:,3))).* min(E_MappedData(:,3)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10); hold on
plot3(ones(size(E_MappedData(:,1))).* min(E_MappedData(:,1)),E_MappedData(:,2),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(E_MappedData(:,1),ones(size(E_MappedData(:,2))).* min(E_MappedData(:,2)),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(E_MappedData(1001,1),E_MappedData(1001,2),E_MappedData(1001,3),100,[0 0.5 0],'diamond','LineWidth',10);
scatter3(E_MappedData(1001,1),E_MappedData(1001,2),min(E_MappedData(:,3)),75,'k','diamond','LineWidth',5);
scatter3(min(E_MappedData(:,1)),E_MappedData(1001,2),E_MappedData(1001,3),75,'k','diamond','LineWidth',5);
scatter3(E_MappedData(1001,1),min(E_MappedData(:,2)),E_MappedData(1001,3),75,'k','diamond','LineWidth',5);
get(gca, 'XTick'); set(gca, 'FontSize', 18); get(gca, 'YTick'); set(gca, 'FontSize', 18);
view(az,el);
text(0.02,0.98,'(c)','Units','Normalized','VerticalAlignment','Top','FontSize',18)

hcb = colorbar('Position',[0.915846994535519,0.381118881118881,0.012021857923497,0.452214452214452],'FontSize',16);
colorTitleHandle = get(hcb,'Title'); titleString = '\eta'; 
set(colorTitleHandle,'String',titleString,'FontSize',30,'Position',[-12.750000000000224,287.9506993006992,0]);
