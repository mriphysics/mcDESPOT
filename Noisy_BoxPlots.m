%%% Generates Figure 7 - fitting of a no exchange model to simulated WML data that has exchange. %%%

T1_S = 1; T1_F = 0.35; T2_S = 0.080; T2_F = 0.015; k_FS = 0:4:20; M0_F = 0.25; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

cmap = colormap(lines(4));
cm = [cmap(1:4,:) ; cmap(1:4,:) ; cmap(1:4,:) ; cmap(1:4,:) ; cmap(1:4,:) ; cmap(1:4,:)];

OnesVector = [ones(1,1000) ; 2*ones(1,1000) ; 3*ones(1,1000) ; 4*ones(1,1000) ; 5*ones(1,1000) ; 6*ones(1,1000) ; 7*ones(1,1000) ; 8*ones(1,1000) ; 9*ones(1,1000) ; 10*ones(1,1000) ; 11*ones(1,1000) ; 12*ones(1,1000) ; 13*ones(1,1000) ; 14*ones(1,1000) ; 15*ones(1,1000) ; 16*ones(1,1000) ; 17*ones(1,1000) ; 18*ones(1,1000) ; 19*ones(1,1000) ; 20*ones(1,1000) ; 21*ones(1,1000) ; 22*ones(1,1000) ; 23*ones(1,1000) ; 24*ones(1,1000)];
Solution_1 = [Solution_B(1,:,1) ; Solution_D(1,:,1) ; Solution_W(1,:,1) ; Solution_Z(1,:,1) ; Solution_B(2,:,1) ; Solution_D(2,:,1) ; Solution_W(2,:,1) ; Solution_Z(2,:,1) ; Solution_B(3,:,1) ; Solution_D(3,:,1) ; Solution_W(3,:,1) ; Solution_Z(3,:,1) ; Solution_B(4,:,1) ; Solution_D(4,:,1) ; Solution_W(4,:,1) ; Solution_Z(4,:,1) ; Solution_B(5,:,1) ; Solution_D(5,:,1) ; Solution_W(5,:,1) ; Solution_Z(5,:,1) ; Solution_B(6,:,1) ; Solution_D(6,:,1) ; Solution_W(6,:,1) ; Solution_Z(6,:,1)].';    
Solution_2 = [Solution_B(1,:,2) ; Solution_D(1,:,2) ; Solution_W(1,:,2) ; Solution_Z(1,:,2) ; Solution_B(2,:,2) ; Solution_D(2,:,2) ; Solution_W(2,:,2) ; Solution_Z(2,:,2) ; Solution_B(3,:,2) ; Solution_D(3,:,2) ; Solution_W(3,:,2) ; Solution_Z(3,:,2) ; Solution_B(4,:,2) ; Solution_D(4,:,2) ; Solution_W(4,:,2) ; Solution_Z(4,:,2) ; Solution_B(5,:,2) ; Solution_D(5,:,2) ; Solution_W(5,:,2) ; Solution_Z(5,:,2) ; Solution_B(6,:,2) ; Solution_D(6,:,2) ; Solution_W(6,:,2) ; Solution_Z(6,:,2)].';
Solution_3 = [Solution_B(1,:,3) ; Solution_D(1,:,3) ; Solution_W(1,:,3) ; Solution_Z(1,:,3) ; Solution_B(2,:,3) ; Solution_D(2,:,3) ; Solution_W(2,:,3) ; Solution_Z(2,:,3) ; Solution_B(3,:,3) ; Solution_D(3,:,3) ; Solution_W(3,:,3) ; Solution_Z(3,:,3) ; Solution_B(4,:,3) ; Solution_D(4,:,3) ; Solution_W(4,:,3) ; Solution_Z(4,:,3) ; Solution_B(5,:,3) ; Solution_D(5,:,3) ; Solution_W(5,:,3) ; Solution_Z(5,:,3) ; Solution_B(6,:,3) ; Solution_D(6,:,3) ; Solution_W(6,:,3) ; Solution_Z(6,:,3)].';
Solution_4 = [Solution_B(1,:,4) ; Solution_D(1,:,4) ; Solution_W(1,:,4) ; Solution_Z(1,:,4) ; Solution_B(2,:,4) ; Solution_D(2,:,4) ; Solution_W(2,:,4) ; Solution_Z(2,:,4) ; Solution_B(3,:,4) ; Solution_D(3,:,4) ; Solution_W(3,:,4) ; Solution_Z(3,:,4) ; Solution_B(4,:,4) ; Solution_D(4,:,4) ; Solution_W(4,:,4) ; Solution_Z(4,:,4) ; Solution_B(5,:,4) ; Solution_D(5,:,4) ; Solution_W(5,:,4) ; Solution_Z(5,:,4) ; Solution_B(6,:,4) ; Solution_D(6,:,4) ; Solution_W(6,:,4) ; Solution_Z(6,:,4)].';
Solution_5 = [Solution_B(1,:,5) ; Solution_D(1,:,5) ; Solution_W(1,:,5) ; Solution_Z(1,:,5) ; Solution_B(2,:,5) ; Solution_D(2,:,5) ; Solution_W(2,:,5) ; Solution_Z(2,:,5) ; Solution_B(3,:,5) ; Solution_D(3,:,5) ; Solution_W(3,:,5) ; Solution_Z(3,:,5) ; Solution_B(4,:,5) ; Solution_D(4,:,5) ; Solution_W(4,:,5) ; Solution_Z(4,:,5) ; Solution_B(5,:,5) ; Solution_D(5,:,5) ; Solution_W(5,:,5) ; Solution_Z(5,:,5) ; Solution_B(6,:,5) ; Solution_D(6,:,5) ; Solution_W(6,:,5) ; Solution_Z(6,:,5)].';

Labels = {'0','0','0','0','4','4','4','4','8','8','8','8','12','12','12','12','16','16','16','16','20','20','20','20'};

figure(1); subplot(2,1,1)
area([0 4.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k'); hold on
area([8.5 12.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
area([16.5 20.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
BP1 = boxplot(Solution_1,OnesVector,'Colors',cm,'Symbol','','Labels',Labels);
for ih=1:7
    set(BP1(ih,:),'LineWidth',2);
end
ylabel('T_{1S} (s)','FontSize',16); hline(T1_S,'k--'); ylim([0.75 1.4])
xticks([]); xticklabels({})
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20); grid on; grid minor;
hLegend = legend(findall(gca,'Tag','Box'), {'B1','B2','B3','B4'},'Position',[0.732422588410056,0.768849862506305,0.169945353302148,0.041379309276055]); hLegend.FontSize = 20; hLegend.Orientation = 'horizontal';
box_vars = findall(gca,'Tag','Box');
hLegend = legend(box_vars([4,3,2,1]), {'B1','B2','B3','B4'}); legend('boxoff')
text(1.75,1.3,'k_{FS} = 0s^{-1}','FontSize',20)
text(5.75,1.3,'k_{FS} = 4s^{-1}','FontSize',20)
text(9.75,1.3,'k_{FS} = 8s^{-1}','FontSize',20)
text(13.75,1.3,'k_{FS} = 12s^{-1}','FontSize',20)
text(17.75,1.3,'k_{FS} = 16s^{-1}','FontSize',20)
text(21.75,1.3,'k_{FS} = 20s^{-1}','FontSize',20)
title('Monte Carlo Simulations - Data Including Exchange and Model Excluding Exchange','FontSize',22,'FontWeight','Bold');

subplot(2,1,2)
area([0 4.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k'); hold on
area([8.5 12.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
area([16.5 20.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
BP2 = boxplot(Solution_2,OnesVector,'Colors',cm,'Symbol','','Labels',Labels); 
for ih=1:7
    set(BP2(ih,:),'LineWidth',2);
end
ylabel('T_{1F} (s)','FontSize',20); hline(T1_F,'k--');
xticks([]); xticklabels({})
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20); grid on; grid minor;
text(1.75,0.65,'k_{FS} = 0s^{-1}','FontSize',20)
text(5.75,0.65,'k_{FS} = 4s^{-1}','FontSize',20)
text(9.75,0.65,'k_{FS} = 8s^{-1}','FontSize',20)
text(13.75,0.65,'k_{FS} = 12s^{-1}','FontSize',20)
text(17.75,0.65,'k_{FS} = 16s^{-1}','FontSize',20)
text(21.75,0.65,'k_{FS} = 20s^{-1}','FontSize',20)

figure(2); subplot(2,1,1)
area([0 4.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k'); hold on
area([8.5 12.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
area([16.5 20.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
BP3 = boxplot(Solution_3,OnesVector,'Colors',cm,'Symbol','','Labels',Labels); 
for ih=1:7
    set(BP3(ih,:),'LineWidth',2);
end
ylabel('T_{2S} (s)','FontSize',20); hline(T2_S,'k--');
xticks([]); xticklabels({})
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20); grid on; grid minor;
text(1.75,0.065,'k_{FS} = 0s^{-1}','FontSize',20)
text(5.75,0.065,'k_{FS} = 4s^{-1}','FontSize',20)
text(9.75,0.065,'k_{FS} = 8s^{-1}','FontSize',20)
text(13.75,0.065,'k_{FS} = 12s^{-1}','FontSize',20)
text(17.75,0.065,'k_{FS} = 16s^{-1}','FontSize',20)
text(21.75,0.065,'k_{FS} = 20s^{-1}','FontSize',20)

subplot(2,1,2)
area([0 4.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k'); hold on
area([8.5 12.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
area([16.5 20.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
BP4 = boxplot(Solution_4,OnesVector,'Colors',cm,'Symbol','','Labels',Labels); 
xticks([]); xticklabels({})
for ih=1:7
    set(BP4(ih,:),'LineWidth',2);
end
ylabel('T_{2F} (s)','FontSize',20); hline(T2_F,'k--'); ylim([0 0.05])
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20); grid on; grid minor;
text(1.75,0.045,'k_{FS} = 0s^{-1}','FontSize',20)
text(5.75,0.045,'k_{FS} = 4s^{-1}','FontSize',20)
text(9.75,0.045,'k_{FS} = 8s^{-1}','FontSize',20)
text(13.75,0.045,'k_{FS} = 12s^{-1}','FontSize',20)
text(17.75,0.045,'k_{FS} = 16s^{-1}','FontSize',20)
text(21.75,0.045,'k_{FS} = 20s^{-1}','FontSize',20)

figure(3); subplot(2,1,1)
area([0 4.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k'); hold on
area([8.5 12.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
area([16.5 20.5],[2 2],'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','k','EdgeColor','k');
BP5 = boxplot(Solution_5,OnesVector,'Colors',cm,'Symbol','','Labels',Labels); 
for ih=1:7
    set(BP5(ih,:),'LineWidth',2);
end
ylabel('MWF','FontSize',20); hline(M0_F,'k--');
xticks([]); xticklabels({})
get(gca, 'XTick'); set(gca, 'FontSize', 20); get(gca, 'YTick'); set(gca, 'FontSize', 20); grid on; grid minor;
text(1.75,0.45,'k_{FS} = 0s^{-1}','FontSize',20)
text(5.75,0.45,'k_{FS} = 4s^{-1}','FontSize',20)
text(9.75,0.45,'k_{FS} = 8s^{-1}','FontSize',20)
text(13.75,0.45,'k_{FS} = 12s^{-1}','FontSize',20)
text(17.75,0.45,'k_{FS} = 16s^{-1}','FontSize',20)
text(21.75,0.45,'k_{FS} = 20s^{-1}','FontSize',20)