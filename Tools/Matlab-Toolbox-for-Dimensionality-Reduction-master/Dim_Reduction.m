% Ground-truth.
T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15;
T1 = 1; T2 = 0.1; M0 = 1;

% Modify data to include ground-truth.
E_TSNE_Matrix_Mod = [E_TSNE_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, 0]]; E_Values_Top_Mod = [E_Values_Top ; 1];
NE_TSNE_Matrix_Mod = [NE_TSNE_Matrix ; [T1_F, T1_S, T2_F, T2_S, M0_F, 0]]; NE_Values_Top_Mod = [NE_Values_Top ; 1];
SP_TSNE_Matrix_Mod = [SP_TSNE_Matrix ; [T1, T2, M0]]; SP_Values_Top_Mod = [SP_Values_Top ; 1];

% Perform dimensionality reduction.
[E_MappedData, E_Mapping] = compute_mapping(E_TSNE_Matrix_Mod(:,1:5), 'KernelPCA', 3);
[NE_MappedData, NE_Mapping] = compute_mapping(NE_TSNE_Matrix_Mod(:,1:5), 'KernelPCA', 3);

% Plot results for model comparison.

az = -159.8; el = 25.2;
Expand1 = 0.005; Expand2 = 0.01; Expand3 = 0.01;

subplot(1,3,3)
colormap(magma); plot3k(E_MappedData,'ColorData',E_Values_Top_Mod,'Marker',{'o',5},'ColorBar',false,'Labels',{'','Dim 1','Dim 2','Dim 3',''}); tt = title('2-Pool (Non-Zero Exchange)'); tt.FontSize = 14; axis square;
hold on; scatter3(E_MappedData(1001,1),E_MappedData(1001,2),E_MappedData(1001,3),20,'g','diamond','LineWidth',5);
view(az,el); 

plot3(E_MappedData(:,1),E_MappedData(:,2),ones(size(E_MappedData(:,3))).* (min(E_MappedData(:,3))-Expand3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(E_MappedData(:,1))).* (max(E_MappedData(:,1))+Expand3),E_MappedData(:,2),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(E_MappedData(:,1),ones(size(E_MappedData(:,2))).* (min(E_MappedData(:,2))-Expand3),E_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(E_MappedData(1001,1),E_MappedData(1001,2),min(E_MappedData(:,3))-Expand3,20,'g','diamond','LineWidth',5);
scatter3(max(E_MappedData(:,1))+Expand3,E_MappedData(1001,2),E_MappedData(1001,3),20,'g','diamond','LineWidth',5);
scatter3(E_MappedData(1001,1),min(E_MappedData(:,2))-Expand3,E_MappedData(1001,3),20,'g','diamond','LineWidth',5);

xlim([min(E_MappedData(:,1))-Expand3 max(E_MappedData(:,1))+Expand3]); ylim([min(E_MappedData(:,2))-Expand3 max(E_MappedData(:,2))+Expand3]); zlim([min(E_MappedData(:,3))-Expand3 max(E_MappedData(:,3))+Expand3]); 

subplot(1,3,2)
colormap(magma); plot3k(NE_MappedData,'ColorData',NE_Values_Top_Mod,'Marker',{'o',5},'ColorBar',false,'Labels',{'','Dim 1','Dim 2','Dim 3',''}); tt = title('2-Pool (Zero Exchange)'); tt.FontSize = 14; axis square; 
hold on; scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),NE_MappedData(1001,3),20,'g','diamond','LineWidth',5)
view(az,el)

plot3(NE_MappedData(:,1),NE_MappedData(:,2),ones(size(NE_MappedData(:,3))).* (min(NE_MappedData(:,3)-Expand2)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(NE_MappedData(:,1))).* (max(NE_MappedData(:,1))+Expand2),NE_MappedData(:,2),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(NE_MappedData(:,1),ones(size(NE_MappedData(:,2))).* (min(NE_MappedData(:,2))-Expand2),NE_MappedData(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(NE_MappedData(1001,1),NE_MappedData(1001,2),min(NE_MappedData(:,3))-Expand2,20,'g','diamond','LineWidth',5);
scatter3(max(NE_MappedData(:,1))+Expand2,NE_MappedData(1001,2),NE_MappedData(1001,3),20,'g','diamond','LineWidth',5);
scatter3(NE_MappedData(1001,1),min(NE_MappedData(:,2))-Expand2,NE_MappedData(1001,3),20,'g','diamond','LineWidth',5);

xlim([min(NE_MappedData(:,1))-Expand2 max(NE_MappedData(:,1))+Expand2]); ylim([min(NE_MappedData(:,2))-Expand2 max(NE_MappedData(:,2))+Expand2]); zlim([min(NE_MappedData(:,3))-Expand2 max(NE_MappedData(:,3))+Expand2]); 

subplot(1,3,1)
colormap(magma); plot3k(SP_TSNE_Matrix_Mod,'ColorData',SP_Values_Top_Mod,'Marker',{'o',5},'ColorBar',false,'Labels',{'','T_{1} (s)','T_{2} (s)','M_{0}',''}); tt = title('Single-Pool'); tt.FontSize = 14; axis square;
hold on; scatter3(SP_TSNE_Matrix_Mod(1001,1),SP_TSNE_Matrix_Mod(1001,2),SP_TSNE_Matrix_Mod(1001,3),20,'g','diamond','LineWidth',5)
view(az,el);

plot3(SP_TSNE_Matrix_Mod(:,1),SP_TSNE_Matrix_Mod(:,2),ones(size(SP_TSNE_Matrix_Mod(:,3))).* (min(SP_TSNE_Matrix_Mod(:,3)-Expand1)),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(ones(size(SP_TSNE_Matrix_Mod(:,1))).* (max(SP_TSNE_Matrix_Mod(:,1))+Expand1),SP_TSNE_Matrix_Mod(:,2),SP_TSNE_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
plot3(SP_TSNE_Matrix_Mod(:,1),ones(size(SP_TSNE_Matrix_Mod(:,2))).* (min(SP_TSNE_Matrix_Mod(:,2))-Expand1),SP_TSNE_Matrix_Mod(:,3),'.','Color',[0.8 0.8 0.8],'MarkerSize',10);
scatter3(SP_TSNE_Matrix_Mod(1001,1),SP_TSNE_Matrix_Mod(1001,2),min(SP_TSNE_Matrix_Mod(:,3))-Expand1,20,'g','diamond','LineWidth',5);
scatter3(max(SP_TSNE_Matrix_Mod(:,1))+Expand1,SP_TSNE_Matrix_Mod(1001,2),SP_TSNE_Matrix_Mod(1001,3),20,'g','diamond','LineWidth',5);
scatter3(SP_TSNE_Matrix_Mod(1001,1),min(SP_TSNE_Matrix_Mod(:,2))-Expand1,SP_TSNE_Matrix_Mod(1001,3),20,'g','diamond','LineWidth',5);

xlim([min(SP_TSNE_Matrix_Mod(:,1))-Expand1 max(SP_TSNE_Matrix_Mod(:,1))+Expand1]); ylim([min(SP_TSNE_Matrix_Mod(:,2))-Expand1 max(SP_TSNE_Matrix_Mod(:,2))+Expand1]); zlim([min(SP_TSNE_Matrix_Mod(:,3))-Expand1 max(SP_TSNE_Matrix_Mod(:,3))+Expand1]); 

%% PCA test.

[E_Coeff, E_Score, E_Latent] = pca(E_TSNE_Matrix(:,1:5),'NumComponents',3);
E_Ratio = E_Latent(1)/E_Latent(2);

[NE_Coeff, NE_Score, NE_Latent] = pca(NE_TSNE_Matrix(:,1:5),'NumComponents',3);
NE_Ratio = NE_Latent(1)/NE_Latent(2);

Difference = E_Ratio/NE_Ratio;