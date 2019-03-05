%%% Generates data for old Supplementary Figure 1. Individual parameters are changed using HB as a basis and MC performed. %%%

close all; clear all

% Diff M0F.
T1F_1 = 0.45; T1S_1 = 1.4; T2F_1 = 0.015; T2S_1 = 0.09; M0F_1 = 0.15; kFS_1 = 8;
T1F_2 = 0.45; T1S_2 = 1.4; T2F_2 = 0.015; T2S_2 = 0.09; M0F_2 = M0F_1*0.75; kFS_2 = 8; 
T1F_3 = 0.45; T1S_3 = 1.4; T2F_3 = 0.015; T2S_3 = 0.09; M0F_3 = M0F_1*1.25; kFS_3 = 8;
% Diff kFS. 
T1F_4 = 0.45; T1S_4 = 1.4; T2F_4 = 0.015; T2S_4 = 0.09; M0F_4 = 0.15; kFS_4 = kFS_1*0.75;
T1F_5 = 0.45; T1S_5 = 1.4; T2F_5 = 0.015; T2S_5 = 0.09; M0F_5 = 0.15; kFS_5 = kFS_1*1.25;
% Diff T1S.
T1F_6 = 0.45; T1S_6 = T1S_1*0.75; T2F_6 = 0.015; T2S_6 = 0.09; M0F_6 = 0.15; kFS_6 = 8;
T1F_7 = 0.45; T1S_7 = T1S_1*1.25; T2F_7 = 0.015; T2S_7 = 0.09; M0F_7 = 0.15; kFS_7 = 8;
% Diff T2S.
T1F_8 = 0.45; T1S_8 = 1.4; T2F_8 = 0.015; T2S_8 = T2S_1*0.75; M0F_8 = 0.15; kFS_8 = 8;
T1F_9 = 0.45; T1S_9 = 1.4; T2F_9 = 0.015; T2S_9 = T2S_1*1.25; M0F_9 = 0.15; kFS_9 = 8;
% Diff T1F.
T1F_10 = T1F_1*0.75; T1S_10 = 1.4; T2F_10 = 0.015; T2S_10 = 0.09; M0F_10 = 0.15; kFS_10 = 8;
T1F_11 = T1F_1*1.25; T1S_11 = 1.4; T2F_11 = 0.015; T2S_11 = 0.09; M0F_11 = 0.15; kFS_11 = 8;
%Diff T2F.
T1F_12 = 0.45; T1S_12 = 1.4; T2F_12 = T2F_1*0.75; T2S_12 = 0.09; M0F_12 = 0.15; kFS_12 = 8;
T1F_13 = 0.45; T1S_13 = 1.4; T2F_13 = T2F_1*1.25; T2S_13 = 0.09; M0F_13 = 0.15; kFS_13 = 8;

Solution_1 = MultipleSRC(T1F_1, T1S_1, T2F_1, T2S_1, M0F_1, kFS_1);
disp('Tissue Set 1 Done')
Solution_2 = MultipleSRC(T1F_2, T1S_2, T2F_2, T2S_2, M0F_2, kFS_2);
disp('Tissue Set 2 Done')
Solution_3 = MultipleSRC(T1F_3, T1S_3, T2F_3, T2S_3, M0F_3, kFS_3);
disp('Tissue Set 3 Done')
Solution_4 = MultipleSRC(T1F_4, T1S_4, T2F_4, T2S_4, M0F_4, kFS_4);
disp('Tissue Set 4 Done')
Solution_5 = MultipleSRC(T1F_5, T1S_5, T2F_5, T2S_5, M0F_5, kFS_5);
disp('Tissue Set 5 Done')
Solution_6 = MultipleSRC(T1F_6, T1S_6, T2F_6, T2S_6, M0F_6, kFS_6);
disp('Tissue Set 6 Done')
Solution_7 = MultipleSRC(T1F_7, T1S_7, T2F_7, T2S_7, M0F_7, kFS_7);
disp('Tissue Set 7 Done')
Solution_8 = MultipleSRC(T1F_8, T1S_8, T2F_8, T2S_8, M0F_8, kFS_8);
disp('Tissue Set 8 Done')
Solution_9 = MultipleSRC(T1F_9, T1S_9, T2F_9, T2S_9, M0F_9, kFS_9);
disp('Tissue Set 9 Done')
Solution_10 = MultipleSRC(T1F_10, T1S_10, T2F_10, T2S_10, M0F_10, kFS_10);
disp('Tissue Set 10 Done')
Solution_11 = MultipleSRC(T1F_11, T1S_11, T2F_11, T2S_11, M0F_11, kFS_11);
disp('Tissue Set 11 Done')
Solution_12 = MultipleSRC(T1F_12, T1S_12, T2F_12, T2S_12, M0F_12, kFS_12);
disp('Tissue Set 12 Done')
Solution_13 = MultipleSRC(T1F_13, T1S_13, T2F_13, T2S_13, M0F_13, kFS_13);
disp('Tissue Set 13 Done')

%% MWF estimation comparison.

cm = colormap(lines(3));

subplot(2,3,1);
histogram(Solution_2(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_3(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); 
line([M0F_2 M0F_2],[0 200],'LineStyle','--','Color',cm(1,:),'LineWidth',2) 
line([M0F_1 M0F_1],[0 200],'LineStyle','--','Color',cm(2,:),'LineWidth',2) 
line([M0F_3 M0F_3],[0 200],'LineStyle','--','Color',cm(3,:),'LineWidth',2) 
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
tt = title('Different GT MWF'); tt.FontSize = 18;

subplot(2,3,2);
histogram(Solution_4(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_5(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); line([M0F_1 M0F_1],[0 200],'LineStyle','--','Color',cm(2,:),'LineWidth',2) 
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
tt = title('Different GT k_{FS}'); tt.FontSize = 18;

subplot(2,3,3);
histogram(Solution_6(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_7(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); line([M0F_1 M0F_1],[0 250],'LineStyle','--','Color',cm(2,:),'LineWidth',2)  
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
ll = legend('0.75 \times HBV','HBV','1.25 \times HBV'); ll.FontSize = 18; legend('boxoff')
tt = title('Different GT T_{1S}'); tt.FontSize = 18;

subplot(2,3,4);
histogram(Solution_8(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_9(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); line([M0F_1 M0F_1],[0 200],'LineStyle','--','Color',cm(2,:),'LineWidth',2) 
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
tt = title('Different GT T_{2S}'); tt.FontSize = 18;

subplot(2,3,5);
histogram(Solution_10(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_11(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); line([M0F_1 M0F_1],[0 200],'LineStyle','--','Color',cm(2,:),'LineWidth',2) 
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
tt = title('Different GT T_{1F}'); tt.FontSize = 18;

subplot(2,3,6);
histogram(Solution_12(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(1,:)); hold on;
histogram(Solution_1(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(2,:)); 
histogram(Solution_13(:,5),'BinWidth',0.02,'FaceAlpha',0.7,'EdgeColor','k','FaceColor',cm(3,:));
xlabel('Estimated MWF','FontSize',16); ylabel('Count','FontSize',16); line([M0F_1 M0F_1],[0 200],'LineStyle','--','Color',cm(2,:),'LineWidth',2) 
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor;
tt = title('Different GT T_{2F}'); tt.FontSize = 18;
