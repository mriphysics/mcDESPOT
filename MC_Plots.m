T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 8; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

cm = colormap(cool(6));

%% T1S
figure(1); subplot(3,2,1);
x = [Solution_ASDeoni(:,1);Solution_ASRui(:,1);Solution_Bouhrara(:,1);Solution_Deoni(:,1);Solution_Wood(:,1);Solution_Zhang(:,1);];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('T_{1S} (s)','FontSize',12); hline(T1_S,'k--'); ylim([1 2.5]);
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;

%% T1F
subplot(3,2,2);
x = [Solution_ASDeoni(:,2);Solution_ASRui(:,2);Solution_Bouhrara(:,2);Solution_Deoni(:,2);Solution_Wood(:,2);Solution_Zhang(:,2)];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('T_{1F} (s)','FontSize',12); hline(T1_F,'k--');
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;

%% T2S
subplot(3,2,3);
x = [Solution_ASDeoni(:,3);Solution_ASRui(:,3);Solution_Bouhrara(:,3);Solution_Deoni(:,3);Solution_Wood(:,3);Solution_Zhang(:,3)];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('T_{2S} (s)','FontSize',12); hline(T2_S,'k--');
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;

%% T2F
subplot(3,2,4);
x = [Solution_ASDeoni(:,4);Solution_ASRui(:,4);Solution_Bouhrara(:,4);Solution_Deoni(:,4);Solution_Wood(:,4);Solution_Zhang(:,4)];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('T_{2F} (s)','FontSize',12); hline(T2_F,'k--');
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;

%% M0F
subplot(3,2,5);
x = [Solution_ASDeoni(:,5);Solution_ASRui(:,5);Solution_Bouhrara(:,5);Solution_Deoni(:,5);Solution_Wood(:,5);Solution_Zhang(:,5)];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('M_{0F}','FontSize',12); hline(M0_F,'k--');
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;

%% kFS
subplot(3,2,6);
x = [Solution_ASDeoni(:,6);Solution_ASRui(:,6);Solution_Bouhrara(:,6);Solution_Deoni(:,6);Solution_Wood(:,6);Solution_Zhang(:,6)];
g = [ones(1000,1); 2*ones(1000,1); 3*ones(1000,1); 4*ones(1000,1); 5*ones(1000,1); 6*ones(1000,1)];
BP1 = boxplot(x,g,'Colors',cm,'Symbol','','Labels',{'S1 B1','S2 B1','S3 B1','S3 B2','S3 B3','S3 B4'}); 
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('k_{FS} (s^{-1})','FontSize',12); hline(k_FS,'k--');
get(gca, 'XTick'); set(gca, 'FontSize', 12); get(gca, 'YTick'); set(gca, 'FontSize', 12);
grid on; grid minor;


%% Degeneracy and SRC comparison.

figure(1); subplot(3,1,1); 
cm = colormap(cool(4));

histogram(M0F_Picked_T1, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.6); hold on

[N1,E1] = histcounts(Solution_T1_Bouhrara(:,5)); 
plot(E1(2:end),N1,'o--','Color',cm(1,:),'LineWidth',2)
[N2,E2] = histcounts(Solution_T1_Deoni(:,5));
plot(E2(2:end),N2,'o--','Color',cm(2,:),'LineWidth',2)
[N3,E3] = histcounts(Solution_T1_Wood(:,5)); 
plot(E3(2:end),N3,'o--','Color',cm(3,:),'LineWidth',2)
[N4,E4] = histcounts(Solution_T1_Zhang(:,5)); 
plot(E4(2:end),N4,'o--','Color',cm(4,:),'LineWidth',2)
ll = legend('Top 1000 Solutions','SRC + Set 1','SRC + Set 2','SRC + Set 3','SRC + Set 4');
ll.FontSize = 14; ll.AutoUpdate = 'off'; ll.Location = 'northwest';

line([0.25 0.25],[0 250],'LineStyle','--','Color','k','LineWidth',2)
line([0 0],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2); line([0.5 0.5],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2)
line([0 0],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2)
line([0.001 0.001],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2)
line([1e-7 1e-7],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2); line([0.3 0.3],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2)

tt = title('Tissue 1'); tt.FontSize = 16; grid on; grid minor
xlabel('M_{0F}','FontSize',14); ylabel('Count','FontSize', 14); xlim([0 0.55]);

subplot(3,1,2);
histogram(M0F_Picked_T2, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.6); hold on

[N4,E4] = histcounts(Solution_T2_Bouhrara(:,5)); 
plot(E4(2:end),N4,'o--','Color',cm(1,:),'LineWidth',2)
[N5,E5] = histcounts(Solution_T2_Deoni(:,5));
plot(E5(2:end),N5,'o--','Color',cm(2,:),'LineWidth',2)
[N6,E6] = histcounts(Solution_T2_Wood(:,5)); 
plot(E6(2:end),N6,'o--','Color',cm(3,:),'LineWidth',2)
[N7,E7] = histcounts(Solution_T2_Zhang(:,5)); 
plot(E7(2:end),N7,'o--','Color',cm(4,:),'LineWidth',2)

line([0.175 0.175],[0 250],'LineStyle','--','Color','k','LineWidth',2) 
line([0 0],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2); line([0.5 0.5],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2)
line([0 0],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2)
line([0.001 0.001],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2)
line([1e-7 1e-7],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2); line([0.3 0.3],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2)

tt = title('Tissue 2'); tt.FontSize = 16; grid on; grid minor
xlabel('M_{0F}','FontSize',14); ylabel('Count','FontSize', 14); xlim([0 0.55]);

subplot(3,1,3);
histogram(M0F_Picked_T3, 'EdgeColor', 'k', 'FaceColor', 'k', 'FaceAlpha', 0.6); hold on

[N8,E8] = histcounts(Solution_T3_Bouhrara(:,5)); 
plot(E8(2:end),N8,'o--','Color',cm(1,:),'LineWidth',2)
[N9,E9] = histcounts(Solution_T3_Deoni(:,5));
plot(E9(2:end),N9,'o--','Color',cm(2,:),'LineWidth',2)
[N10,E10] = histcounts(Solution_T3_Wood(:,5)); 
plot(E10(2:end),N10,'o--','Color',cm(3,:),'LineWidth',2)
[N11,E11] = histcounts(Solution_T3_Zhang(:,5)); 
plot(E11(2:end),N11,'o--','Color',cm(4,:),'LineWidth',2)

line([0.1 0.1],[0 250],'LineStyle','--','Color','k','LineWidth',2)
line([0 0],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2); line([0.5 0.5],[0 250],'LineStyle',':','Color',cm(1,:),'LineWidth',2)
line([0 0],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(2,:),'LineWidth',2)
line([0.001 0.001],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2); line([0.35 0.35],[0 250],'LineStyle',':','Color',cm(3,:),'LineWidth',2)
line([1e-7 1e-7],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2); line([0.3 0.3],[0 250],'LineStyle',':','Color',cm(4,:),'LineWidth',2)

tt = title('Tissue 3'); tt.FontSize = 16; grid on; grid minor
xlabel('M_{0F}','FontSize',14); ylabel('Count','FontSize', 14); xlim([0 0.55]);
