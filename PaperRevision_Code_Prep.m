%%% Test stability of SRC estimates for different no. realisations.  %%%

close all; clear all;

load('Fig5_Data.mat')

Bouhrara_Mean = mean(Solution_Bouhrara);
Bouhrara_Mean_US900 = mean(Solution_Bouhrara(1:900,:));
Bouhrara_Mean_US800 = mean(Solution_Bouhrara(1:800,:));
Bouhrara_Mean_US700 = mean(Solution_Bouhrara(1:700,:));
Bouhrara_Mean_US600 = mean(Solution_Bouhrara(1:600,:));
Bouhrara_Mean_US500 = mean(Solution_Bouhrara(1:500,:));
Bouhrara_Mean_US400 = mean(Solution_Bouhrara(1:400,:));
Bouhrara_Mean_US300 = mean(Solution_Bouhrara(1:300,:));
Bouhrara_Mean_US200 = mean(Solution_Bouhrara(1:200,:));
Bouhrara_Mean_US100 = mean(Solution_Bouhrara(1:100,:));

Bouhrara_SD = std(Solution_Bouhrara);
Bouhrara_SD_US900 = std(Solution_Bouhrara(1:900,:));
Bouhrara_SD_US800 = std(Solution_Bouhrara(1:800,:));
Bouhrara_SD_US700 = std(Solution_Bouhrara(1:700,:));
Bouhrara_SD_US600 = std(Solution_Bouhrara(1:600,:));
Bouhrara_SD_US500 = std(Solution_Bouhrara(1:500,:));
Bouhrara_SD_US400 = std(Solution_Bouhrara(1:400,:));
Bouhrara_SD_US300 = std(Solution_Bouhrara(1:300,:));
Bouhrara_SD_US200 = std(Solution_Bouhrara(1:200,:));
Bouhrara_SD_US100 = std(Solution_Bouhrara(1:100,:));

Bouhrara_MeanVector = [Bouhrara_Mean_US100;Bouhrara_Mean_US200;Bouhrara_Mean_US300;Bouhrara_Mean_US400;Bouhrara_Mean_US500;Bouhrara_Mean_US600;Bouhrara_Mean_US700;Bouhrara_Mean_US800;Bouhrara_Mean_US900;Bouhrara_Mean];
Bouhrara_SDVector = [Bouhrara_SD_US100;Bouhrara_SD_US200;Bouhrara_SD_US300;Bouhrara_SD_US400;Bouhrara_SD_US500;Bouhrara_SD_US600;Bouhrara_SD_US700;Bouhrara_SD_US800;Bouhrara_SD_US900;Bouhrara_SD];

Deoni_Mean = mean(Solution_Deoni);
Deoni_Mean_US900 = mean(Solution_Deoni(1:900,:));
Deoni_Mean_US800 = mean(Solution_Deoni(1:800,:));
Deoni_Mean_US700 = mean(Solution_Deoni(1:700,:));
Deoni_Mean_US600 = mean(Solution_Deoni(1:600,:));
Deoni_Mean_US500 = mean(Solution_Deoni(1:500,:));
Deoni_Mean_US400 = mean(Solution_Deoni(1:400,:));
Deoni_Mean_US300 = mean(Solution_Deoni(1:300,:));
Deoni_Mean_US200 = mean(Solution_Deoni(1:200,:));
Deoni_Mean_US100 = mean(Solution_Deoni(1:100,:));

Deoni_SD = std(Solution_Deoni);
Deoni_SD_US900 = std(Solution_Deoni(1:900,:));
Deoni_SD_US800 = std(Solution_Deoni(1:800,:));
Deoni_SD_US700 = std(Solution_Deoni(1:700,:));
Deoni_SD_US600 = std(Solution_Deoni(1:600,:));
Deoni_SD_US500 = std(Solution_Deoni(1:500,:));
Deoni_SD_US400 = std(Solution_Deoni(1:400,:));
Deoni_SD_US300 = std(Solution_Deoni(1:300,:));
Deoni_SD_US200 = std(Solution_Deoni(1:200,:));
Deoni_SD_US100 = std(Solution_Deoni(1:100,:));

Deoni_MeanVector = [Deoni_Mean_US100;Deoni_Mean_US200;Deoni_Mean_US300;Deoni_Mean_US400;Deoni_Mean_US500;Deoni_Mean_US600;Deoni_Mean_US700;Deoni_Mean_US800;Deoni_Mean_US900;Deoni_Mean];
Deoni_SDVector = [Deoni_SD_US100;Deoni_SD_US200;Deoni_SD_US300;Deoni_SD_US400;Deoni_SD_US500;Deoni_SD_US600;Deoni_SD_US700;Deoni_SD_US800;Deoni_SD_US900;Deoni_SD];

Wood_Mean = mean(Solution_Wood);
Wood_Mean_US900 = mean(Solution_Wood(1:900,:));
Wood_Mean_US800 = mean(Solution_Wood(1:800,:));
Wood_Mean_US700 = mean(Solution_Wood(1:700,:));
Wood_Mean_US600 = mean(Solution_Wood(1:600,:));
Wood_Mean_US500 = mean(Solution_Wood(1:500,:));
Wood_Mean_US400 = mean(Solution_Wood(1:400,:));
Wood_Mean_US300 = mean(Solution_Wood(1:300,:));
Wood_Mean_US200 = mean(Solution_Wood(1:200,:));
Wood_Mean_US100 = mean(Solution_Wood(1:100,:));

Wood_SD = std(Solution_Wood);
Wood_SD_US900 = std(Solution_Wood(1:900,:));
Wood_SD_US800 = std(Solution_Wood(1:800,:));
Wood_SD_US700 = std(Solution_Wood(1:700,:));
Wood_SD_US600 = std(Solution_Wood(1:600,:));
Wood_SD_US500 = std(Solution_Wood(1:500,:));
Wood_SD_US400 = std(Solution_Wood(1:400,:));
Wood_SD_US300 = std(Solution_Wood(1:300,:));
Wood_SD_US200 = std(Solution_Wood(1:200,:));
Wood_SD_US100 = std(Solution_Wood(1:100,:));

Wood_MeanVector = [Wood_Mean_US100;Wood_Mean_US200;Wood_Mean_US300;Wood_Mean_US400;Wood_Mean_US500;Wood_Mean_US600;Wood_Mean_US700;Wood_Mean_US800;Wood_Mean_US900;Wood_Mean];
Wood_SDVector = [Wood_SD_US100;Wood_SD_US200;Wood_SD_US300;Wood_SD_US400;Wood_SD_US500;Wood_SD_US600;Wood_SD_US700;Wood_SD_US800;Wood_SD_US900;Wood_SD];

Zhang_Mean = mean(Solution_Zhang);
Zhang_Mean_US900 = mean(Solution_Zhang(1:900,:));
Zhang_Mean_US800 = mean(Solution_Zhang(1:800,:));
Zhang_Mean_US700 = mean(Solution_Zhang(1:700,:));
Zhang_Mean_US600 = mean(Solution_Zhang(1:600,:));
Zhang_Mean_US500 = mean(Solution_Zhang(1:500,:));
Zhang_Mean_US400 = mean(Solution_Zhang(1:400,:));
Zhang_Mean_US300 = mean(Solution_Zhang(1:300,:));
Zhang_Mean_US200 = mean(Solution_Zhang(1:200,:));
Zhang_Mean_US100 = mean(Solution_Zhang(1:100,:));

Zhang_SD = std(Solution_Zhang);
Zhang_SD_US900 = std(Solution_Zhang(1:900,:));
Zhang_SD_US800 = std(Solution_Zhang(1:800,:));
Zhang_SD_US700 = std(Solution_Zhang(1:700,:));
Zhang_SD_US600 = std(Solution_Zhang(1:600,:));
Zhang_SD_US500 = std(Solution_Zhang(1:500,:));
Zhang_SD_US400 = std(Solution_Zhang(1:400,:));
Zhang_SD_US300 = std(Solution_Zhang(1:300,:));
Zhang_SD_US200 = std(Solution_Zhang(1:200,:));
Zhang_SD_US100 = std(Solution_Zhang(1:100,:));

Zhang_MeanVector = [Zhang_Mean_US100;Zhang_Mean_US200;Zhang_Mean_US300;Zhang_Mean_US400;Zhang_Mean_US500;Zhang_Mean_US600;Zhang_Mean_US700;Zhang_Mean_US800;Zhang_Mean_US900;Zhang_Mean];
Zhang_SDVector = [Zhang_SD_US100;Zhang_SD_US200;Zhang_SD_US300;Zhang_SD_US400;Zhang_SD_US500;Zhang_SD_US600;Zhang_SD_US700;Zhang_SD_US800;Zhang_SD_US900;Zhang_SD];

ASDeoni_Mean = mean(Solution_ASDeoni);
ASDeoni_Mean_US900 = mean(Solution_ASDeoni(1:900,:));
ASDeoni_Mean_US800 = mean(Solution_ASDeoni(1:800,:));
ASDeoni_Mean_US700 = mean(Solution_ASDeoni(1:700,:));
ASDeoni_Mean_US600 = mean(Solution_ASDeoni(1:600,:));
ASDeoni_Mean_US500 = mean(Solution_ASDeoni(1:500,:));
ASDeoni_Mean_US400 = mean(Solution_ASDeoni(1:400,:));
ASDeoni_Mean_US300 = mean(Solution_ASDeoni(1:300,:));
ASDeoni_Mean_US200 = mean(Solution_ASDeoni(1:200,:));
ASDeoni_Mean_US100 = mean(Solution_ASDeoni(1:100,:));

ASDeoni_SD = std(Solution_ASDeoni);
ASDeoni_SD_US900 = std(Solution_ASDeoni(1:900,:));
ASDeoni_SD_US800 = std(Solution_ASDeoni(1:800,:));
ASDeoni_SD_US700 = std(Solution_ASDeoni(1:700,:));
ASDeoni_SD_US600 = std(Solution_ASDeoni(1:600,:));
ASDeoni_SD_US500 = std(Solution_ASDeoni(1:500,:));
ASDeoni_SD_US400 = std(Solution_ASDeoni(1:400,:));
ASDeoni_SD_US300 = std(Solution_ASDeoni(1:300,:));
ASDeoni_SD_US200 = std(Solution_ASDeoni(1:200,:));
ASDeoni_SD_US100 = std(Solution_ASDeoni(1:100,:));

ASDeoni_MeanVector = [ASDeoni_Mean_US100;ASDeoni_Mean_US200;ASDeoni_Mean_US300;ASDeoni_Mean_US400;ASDeoni_Mean_US500;ASDeoni_Mean_US600;ASDeoni_Mean_US700;ASDeoni_Mean_US800;ASDeoni_Mean_US900;ASDeoni_Mean];
ASDeoni_SDVector = [ASDeoni_SD_US100;ASDeoni_SD_US200;ASDeoni_SD_US300;ASDeoni_SD_US400;ASDeoni_SD_US500;ASDeoni_SD_US600;ASDeoni_SD_US700;ASDeoni_SD_US800;ASDeoni_SD_US900;ASDeoni_SD];

ASRui_Mean = mean(Solution_ASRui);
ASRui_Mean_US900 = mean(Solution_ASRui(1:900,:));
ASRui_Mean_US800 = mean(Solution_ASRui(1:800,:));
ASRui_Mean_US700 = mean(Solution_ASRui(1:700,:));
ASRui_Mean_US600 = mean(Solution_ASRui(1:600,:));
ASRui_Mean_US500 = mean(Solution_ASRui(1:500,:));
ASRui_Mean_US400 = mean(Solution_ASRui(1:400,:));
ASRui_Mean_US300 = mean(Solution_ASRui(1:300,:));
ASRui_Mean_US200 = mean(Solution_ASRui(1:200,:));
ASRui_Mean_US100 = mean(Solution_ASRui(1:100,:));

ASRui_SD = std(Solution_ASRui);
ASRui_SD_US900 = std(Solution_ASRui(1:900,:));
ASRui_SD_US800 = std(Solution_ASRui(1:800,:));
ASRui_SD_US700 = std(Solution_ASRui(1:700,:));
ASRui_SD_US600 = std(Solution_ASRui(1:600,:));
ASRui_SD_US500 = std(Solution_ASRui(1:500,:));
ASRui_SD_US400 = std(Solution_ASRui(1:400,:));
ASRui_SD_US300 = std(Solution_ASRui(1:300,:));
ASRui_SD_US200 = std(Solution_ASRui(1:200,:));
ASRui_SD_US100 = std(Solution_ASRui(1:100,:));

ASRui_MeanVector = [ASRui_Mean_US100;ASRui_Mean_US200;ASRui_Mean_US300;ASRui_Mean_US400;ASRui_Mean_US500;ASRui_Mean_US600;ASRui_Mean_US700;ASRui_Mean_US800;ASRui_Mean_US900;ASRui_Mean];
ASRui_SDVector = [ASRui_SD_US100;ASRui_SD_US200;ASRui_SD_US300;ASRui_SD_US400;ASRui_SD_US500;ASRui_SD_US600;ASRui_SD_US700;ASRui_SD_US800;ASRui_SD_US900;ASRui_SD];

USFactors = 100:100:1000;

figure(1); subplot(3,2,1); 
plot(USFactors,Bouhrara_MeanVector(:,1),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,1),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,1),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,1),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,1),'-o','LineWidth',2);
plot(USFactors,ASRui_MeanVector(:,1),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('Mean T_{1S} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])
ll = legend('S1+B1','S1+B2','S1+B3','S1+B4','S2+B1','S3+B1'); ll.FontSize = 16; legend 'boxoff'; ll.Orientation = 'horizontal';

subplot(3,2,2); 
plot(USFactors,Bouhrara_MeanVector(:,2),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,2),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,2),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,2),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,2),'-o','LineWidth',2);
plot(USFactors,ASRui_MeanVector(:,2),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('Mean T_{1F} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,3); 
plot(USFactors,Bouhrara_MeanVector(:,3),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,3),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,3),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,3),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,3),'-o','LineWidth',2);
plot(USFactors,ASRui_MeanVector(:,3),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('Mean T_{2S} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,4); 
plot(USFactors,Bouhrara_MeanVector(:,4),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,4),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,4),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,4),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,4),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('Mean T_{2F} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,5); 
plot(USFactors,Bouhrara_MeanVector(:,5),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,5),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,5),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,5),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,5),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16);ylabel('Mean MWF','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,6); 
plot(USFactors,Bouhrara_MeanVector(:,6),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_MeanVector(:,6),'-o','LineWidth',2);
plot(USFactors,Wood_MeanVector(:,6),'-o','LineWidth',2);
plot(USFactors,Zhang_MeanVector(:,6),'-o','LineWidth',2);
plot(USFactors,ASDeoni_MeanVector(:,6),'-o','LineWidth',2);
plot(USFactors,ASRui_MeanVector(:,6),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('Mean k_{FS} (s^{-1})','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

figure(2); subplot(3,2,1); 
plot(USFactors,Bouhrara_SDVector(:,1),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,1),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,1),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,1),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,1),'-o','LineWidth',2);
plot(USFactors,ASRui_SDVector(:,1),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('SD T_{1S} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])
ll = legend('S1+B1','S1+B2','S1+B3','S1+B4','S2+B1','S3+B1'); ll.FontSize = 16; legend 'boxoff'; ll.Orientation = 'horizontal';

subplot(3,2,2); 
plot(USFactors,Bouhrara_SDVector(:,2),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,2),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,2),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,2),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,2),'-o','LineWidth',2);
plot(USFactors,ASRui_SDVector(:,2),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('SD T_{1F} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,3); 
plot(USFactors,Bouhrara_SDVector(:,3),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,3),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,3),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,3),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,3),'-o','LineWidth',2);
plot(USFactors,ASRui_SDVector(:,3),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('SD T_{2S} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,4); 
plot(USFactors,Bouhrara_SDVector(:,4),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,4),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,4),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,4),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,4),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('SD T_{2F} (s)','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,5); 
plot(USFactors,Bouhrara_SDVector(:,5),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,5),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,5),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,5),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,5),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16);ylabel('SD MWF','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])

subplot(3,2,6); 
plot(USFactors,Bouhrara_SDVector(:,6),'-o','LineWidth',2); hold on;
plot(USFactors,Deoni_SDVector(:,6),'-o','LineWidth',2);
plot(USFactors,Wood_SDVector(:,6),'-o','LineWidth',2);
plot(USFactors,Zhang_SDVector(:,6),'-o','LineWidth',2);
plot(USFactors,ASDeoni_SDVector(:,6),'-o','LineWidth',2);
plot(USFactors,ASRui_SDVector(:,6),'-o','LineWidth',2);
xlabel('No. Realisations','FontSize',16); ylabel('SD k_{FS} (s^{-1})','FontSize',16);
get(gca, 'XTick'); set(gca, 'FontSize', 16); get(gca, 'YTick'); set(gca, 'FontSize', 16);
grid on; grid minor; xlim([100 1000])