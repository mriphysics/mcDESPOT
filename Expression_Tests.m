%% Investigation into MWF 1.

close all
clear all

T1_S = 1.15; T1_F = 0.4; k_FS = 9; T1_B = 1; k_B = 5; k_FB = 5; k_SB = 5;
M0_F = 0:0.01:0.5; 
M0_B = [0.1 0.2 0.3 0.4 0.5]; 

TR = 5.2e-3; dM_zB = 0;
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; B1 = 13e-6;
T_RF = deg2rad(50)/(Gamma * B1); W = (pi/T_RF) * (Gamma * B1)^2 * T_RF * G * (T_RF/TR);

k_FSapp = zeros(length(M0_F),length(M0_B)); k_SFapp = zeros(length(M0_F),length(M0_B));
T1_Fapp = zeros(length(M0_F),length(M0_B)); T1_Sapp = zeros(length(M0_F),length(M0_B));
M_0Fapp = zeros(length(M0_F),length(M0_B)); M_0Sapp = zeros(length(M0_F),length(M0_B));

k_FSapp_NA = zeros(length(M0_F),length(M0_B)); k_SFapp_NA = zeros(length(M0_F),length(M0_B));
T1_Fapp_NA = zeros(length(M0_F),length(M0_B)); T1_Sapp_NA = zeros(length(M0_F),length(M0_B));
M_0Fapp_NA = zeros(length(M0_F),length(M0_B)); M_0Sapp_NA = zeros(length(M0_F),length(M0_B));

MWF_app = zeros(length(M0_F),length(M0_B));
MWF_app_NA = zeros(length(M0_F),length(M0_B));
M0_S = zeros(length(M0_F),length(M0_B)); MWF = zeros(length(M0_F),length(M0_B));

for ss = 1:length(M0_F)
    
    for tt = 1:length(M0_B)
    
    MWF(ss,tt) = M0_F(ss)/(1-M0_B(tt)); % MWF(ss,tt) = M0_F(ss)/(M0_F(ss) + M0_S(ss,tt));
    M0_S(ss,tt) = 1 - (M0_F(ss) + M0_B(tt));
    %k_SF = (M0_F * k_FS)/M0_S;

    k_FSapp(ss,tt) = (M0_S(ss,tt)*k_B^2)/(M0_B(tt)*W + M0_F(ss)*k_B + M0_S(ss,tt)*k_B + M0_B(tt)*R1_B);
    k_SFapp(ss,tt) = (M0_F(ss)*k_B^2)/(M0_B(tt)*W + M0_F(ss)*k_B + M0_S(ss,tt)*k_B + M0_B(tt)*R1_B);
    T1_Fapp(ss,tt) = 1/((M0_B(tt)*R1_B*R1_F + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss,tt)*R1_F*k_B + M0_B(tt)*W*k_B)/(M0_B(tt)*W + M0_F(ss)*k_B + M0_S(ss,tt)*k_B + M0_B(tt)*R1_B));
    T1_Sapp(ss,tt) = 1/((M0_B(tt)*R1_B*R1_S + M0_B(tt)*R1_S*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss,tt)*R1_S*k_B + M0_B(tt)*W*k_B)/(M0_B(tt)*W + M0_F(ss)*k_B + M0_S(ss,tt)*k_B + M0_B(tt)*R1_B));
    M_0Fapp(ss,tt) = (M0_F(ss)*(M0_B(tt)*R1_B*R1_F - dM_zB*k_B + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss,tt)*R1_F*k_B))/(M0_B(tt)*R1_B*R1_F + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss,tt)*R1_F*k_B + M0_B(tt)*W*k_B);
    M_0Sapp(ss,tt) = (M0_S(ss,tt)*(M0_B(tt)*R1_B*R1_S - dM_zB*k_B + M0_B(tt)*R1_S*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss,tt)*R1_S*k_B))/(M0_B(tt)*R1_B*R1_S + M0_B(tt)*R1_S*W + M0_B(tt)*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss,tt)*R1_S*k_B + M0_B(tt)*W*k_B);
    MWF_app(ss,tt) = M_0Fapp(ss,tt)/(M_0Fapp(ss,tt) + M_0Sapp(ss,tt));
    
    k_SFapp_NA(ss,tt) = (M0_F(ss)*k_FS)/M0_S(ss,tt) + (M0_F(ss)*k_FB*k_SB)/(M0_B(tt)*W + M0_F(ss)*k_FB + M0_S(ss,tt)*k_SB + M0_B(tt)*R1_B);
    k_FSapp_NA(ss,tt) = k_FS + (M0_S(ss,tt)*k_FB*k_SB)/(M0_B(tt)*W + M0_F(ss)*k_FB + M0_S(ss,tt)*k_SB + M0_B(tt)*R1_B);
    T1_Fapp_NA(ss,tt) = 1/((M0_B(tt)*R1_B*R1_F + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_FB + M0_F(ss)*R1_F*k_FB + M0_S(ss,tt)*R1_F*k_SB + M0_B(tt)*W*k_FB)/(M0_B(tt)*W + M0_F(ss)*k_FB + M0_S(ss,tt)*k_SB + M0_B(tt)*R1_B));
    T1_Sapp_NA(ss,tt) = 1/((M0_B(tt)*R1_B*R1_S + M0_B(tt)*R1_S*W + M0_F(ss)*R1_S*k_FB + M0_B(tt)*R1_B*k_SB + M0_S(ss,tt)*R1_S*k_SB + M0_B(tt)*W*k_SB)/(M0_B(tt)*W + M0_F(ss)*k_FB + M0_S(ss,tt)*k_SB + M0_B(tt)*R1_B));
    M_0Fapp_NA(ss,tt) = (M0_F(ss)*(M0_B(tt)*R1_B*R1_F - dM_zB*k_FB + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_FB + M0_F(ss)*R1_F*k_FB + M0_S(ss,tt)*R1_F*k_SB))/(M0_B(tt)*R1_B*R1_F + M0_B(tt)*R1_F*W + M0_B(tt)*R1_B*k_FB + M0_F(ss)*R1_F*k_FB + M0_S(ss,tt)*R1_F*k_SB + M0_B(tt)*W*k_FB);
    M_0Sapp_NA(ss,tt) = (M0_S(ss,tt)*(M0_B(tt)*R1_B*R1_S - dM_zB*k_SB + M0_B(tt)*R1_S*W + M0_F(ss)*R1_S*k_FB + M0_B(tt)*R1_B*k_SB + M0_S(ss,tt)*R1_S*k_SB))/(M0_B(tt)*R1_B*R1_S + M0_B(tt)*R1_S*W + M0_F(ss)*R1_S*k_FB + M0_B(tt)*R1_B*k_SB + M0_S(ss,tt)*R1_S*k_SB + M0_B(tt)*W*k_SB);
    MWF_app_NA(ss,tt) = M_0Fapp_NA(ss,tt)/(M_0Fapp_NA(ss,tt) + M_0Sapp_NA(ss,tt));
        
    end
    
end

% kFS_PC = ((k_FSapp_NA - k_FS)/k_FS)*100;
% kSF_PC = ((k_SFapp_NA - k_SF)/k_SF)*100;
% T1S_PC = ((T1_Sapp_NA - T1_S)/T1_S)*100;
% T1F_PC = ((T1_Fapp_NA - T1_F)/T1_F)*100;
% M0F_PC = ((M_0Fapp_NA - M0_F)/M0_F)*100;
% M0S_PC = ((M_0Sapp_NA - M0_S)/M0_S)*100;
% MWF_PC = ((MWF_app_NA - MWF)/MWF)*100;
% 
% kFS_PD = ((k_FSapp - k_FSapp_NA)/k_FSapp)*100;
% kSF_PD = ((k_SFapp - k_SFapp_NA)/k_SFapp)*100;

figure(1)
hold on
subplot(1,3,1)
cm = colormap(cool(length(M0_B)));

for nn = 1:length(M0_B)
    plot(M0_F, MWF_app(:,nn),'LineWidth', 1.5,'color',cm(nn,:))
    hold on
end

xlabel('M_{0F}', 'FontSize', 18)
ylabel('MWF^{app}', 'FontSize', 18)
xlim([0 0.5]) 
Identity = refline(1,0);
Identity.Color = 'k';
Identity.LineStyle = '--';
l = legend('0.1', '0.2', '0.3', '0.4', '0.5');
title(l, 'M_{0B}','FontSize', 18)
l.FontSize = 18;
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
title('(a)','FontSize',18)

subplot(1,3,2)
for tt = 1:length(M0_B)
    plot(MWF(:,tt), MWF_app(:,tt),'LineWidth', 1.5,'color',cm(tt,:))
    hold on
end
xlabel('MWF', 'FontSize', 18)
ylabel('MWF^{app}', 'FontSize', 18)
xlim([0 0.5]) 
Identity = refline(1,0);
Identity.Color = 'k';
Identity.LineStyle = '--';
l = legend('0.1', '0.2', '0.3', '0.4', '0.5');
title(l, 'M_{0B}','FontSize', 18)
l.FontSize = 18;
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
title('(b)','FontSize',18)

%% Investigation into MWF 2.

clear all

T1_S = 1.15; T1_F = 0.4; T1_B = 1; M0_B = 0.2; k_B = 5;
M0_F = 0:0.01:0.5; 
TR_SSFP = 5.2e-3; dM_zB = 0; %TR = 5.2e-3; 
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; 
B1rms = 0.5e-6:0.5e-6:3e-6;  

k_FSapp = zeros(length(B1rms),length(M0_F)); k_SFapp = zeros(length(B1rms),length(M0_F));
T1_Fapp = zeros(length(B1rms),length(M0_F)); T1_Sapp = zeros(length(B1rms),length(M0_F));
M_0Fapp = zeros(length(B1rms),length(M0_F)); M_0Sapp = zeros(length(B1rms),length(M0_F));
%T_RF = zeros(length(B1),1); B1rms = zeros(length(B1),1); W = zeros(length(B1),1);
M0_S = zeros(length(M0_F),1); MWF = zeros(length(M0_F),1);


for ss = 1:length(M0_F)
   
    M0_S(ss) = 1 - (M0_F(ss) + M0_B);
    MWF(ss) = M0_F(ss)/(1-M0_B);

    for tt = 1:length(B1rms)
    
        %T_RF(tt) = deg2rad(50)/(Gamma * B1(tt));
        %B1rms(tt) = sqrt((B1(tt)^2*T_RF(tt))/TR);
        W(tt) = pi * G * (Gamma * B1rms(tt))^2;

        k_FSapp(ss,tt) = (M0_S(ss)*k_B^2)/(M0_B*W(tt) + M0_F(ss)*k_B + M0_S(ss)*k_B + M0_B*R1_B);
        k_SFapp(ss,tt) = (M0_F(ss)*k_B^2)/(M0_B*W(tt) + M0_F(ss)*k_B + M0_S(ss)*k_B + M0_B*R1_B);
        T1_Fapp(ss,tt) = 1/((M0_B*R1_B*R1_F + M0_B*R1_F*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss)*R1_F*k_B + M0_B*W(tt)*k_B)/(M0_B*W(tt) + M0_F(ss)*k_B + M0_S(ss)*k_B + M0_B*R1_B));
        T1_Sapp(ss,tt) = 1/((M0_B*R1_B*R1_S + M0_B*R1_S*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss)*R1_S*k_B + M0_B*W(tt)*k_B)/(M0_B*W(tt) + M0_F(ss)*k_B + M0_S(ss)*k_B + M0_B*R1_B));
        M_0Fapp(ss,tt) = (M0_F(ss)*(M0_B*R1_B*R1_F - dM_zB*k_B + M0_B*R1_F*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss)*R1_F*k_B))/(M0_B*R1_B*R1_F + M0_B*R1_F*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_F*k_B + M0_S(ss)*R1_F*k_B + M0_B*W(tt)*k_B);
        M_0Sapp(ss,tt) = (M0_S(ss)*(M0_B*R1_B*R1_S - dM_zB*k_B + M0_B*R1_S*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss)*R1_S*k_B))/(M0_B*R1_B*R1_S + M0_B*R1_S*W(tt) + M0_B*R1_B*k_B + M0_F(ss)*R1_S*k_B + M0_S(ss)*R1_S*k_B + M0_B*W(tt)*k_B);
        MWF_app(ss,tt) = M_0Fapp(ss,tt)/(M_0Fapp(ss,tt) + M_0Sapp(ss,tt));
    
    end
    
end

figure(1); subplot(1,3,3)
cm = colormap(cool(length(B1rms)));
for nn = 1:length(B1rms)
    plot(MWF, MWF_app(:,nn),'LineWidth', 1.5,'color',cm(nn,:))
    hold on
end
xlabel('MWF', 'FontSize', 18)
ylabel('MWF^{app}', 'FontSize', 18)
xlim([0 0.6250])
Identity = refline(1,0);
Identity.Color = 'k';
Identity.LineStyle = '--';
l = legend('0.5', '1.0', '1.5', '2.0', '2.5', '3.0');
title(l, 'B_{1,rms} [\muT]','FontSize', 18)
l.FontSize = 18;
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
title('(c)','FontSize',18)

%% Investigation into relaxation times.

clear all

T1_S = 1.15; T1_F = 0.4; T1_B = 1; M0_F = 0.25; M0_B = 0.2; M0_S = 1 - (M0_F + M0_B); k_B = 5;

TR = 5.2e-3; dM_zB = 0;
R1_F = 1/T1_F; R1_B = 1/T1_B; R1_S = 1/T1_S;
G = 1.4e-5; Gamma = 2 * pi * 42.57747892e6; 
B1 = 1e-6:0.5e-6:15e-6;  

k_FSapp = zeros(length(B1),1); k_SFapp = zeros(length(B1),1);
T1_Fapp = zeros(length(B1),1); T1_Sapp = zeros(length(B1),1);
M_0Fapp = zeros(length(B1),1); M_0Sapp = zeros(length(B1),1);
T_RF = zeros(length(B1),1); B1rms = zeros(length(B1),1); W = zeros(length(B1),1);

for ss = 1:length(B1)
    
    T_RF(ss) = deg2rad(50)/(Gamma * B1(ss));
    B1rms(ss) = sqrt((B1(ss)^2*T_RF(ss))/TR);
    W(ss) = pi * G * (Gamma * B1rms(ss))^2;
    
    for tt = 1:length(k_B)
        
        k_FSapp(ss) = (M0_S*k_B^2)/(M0_B*W(ss) + M0_F*k_B + M0_S*k_B + M0_B*R1_B);
        k_SFapp(ss) = (M0_F*k_B^2)/(M0_B*W(ss) + M0_F*k_B + M0_S*k_B + M0_B*R1_B);
        T1_Fapp(ss) = 1/((M0_B*R1_B*R1_F + M0_B*R1_F*W(ss) + M0_B*R1_B*k_B + M0_F*R1_F*k_B + M0_S*R1_F*k_B + M0_B*W(ss)*k_B)/(M0_B*W(ss) + M0_F*k_B + M0_S*k_B + M0_B*R1_B));
        T1_Sapp(ss) = 1/((M0_B*R1_B*R1_S + M0_B*R1_S*W(ss) + M0_B*R1_B*k_B + M0_F*R1_S*k_B + M0_S*R1_S*k_B + M0_B*W(ss)*k_B)/(M0_B*W(ss) + M0_F*k_B + M0_S*k_B + M0_B*R1_B));
        M_0Fapp(ss) = (M0_F*(M0_B*R1_B*R1_F - dM_zB*k_B + M0_B*R1_F*W(ss) + M0_B*R1_B*k_B + M0_F*R1_F*k_B + M0_S*R1_F*k_B))/(M0_B*R1_B*R1_F + M0_B*R1_F*W(ss) + M0_B*R1_B*k_B + M0_F*R1_F*k_B + M0_S*R1_F*k_B + M0_B*W(ss)*k_B);
        M_0Sapp(ss) = (M0_S*(M0_B*R1_B*R1_S - dM_zB*k_B + M0_B*R1_S*W(ss) + M0_B*R1_B*k_B + M0_F*R1_S*k_B + M0_S*R1_S*k_B))/(M0_B*R1_B*R1_S + M0_B*R1_S*W(ss) + M0_B*R1_B*k_B + M0_F*R1_S*k_B + M0_S*R1_S*k_B + M0_B*W(ss)*k_B);
        MWF_app(ss) = M_0Fapp(ss)/(M_0Fapp(ss) + M_0Sapp(ss));
        
    end
    
end

% Plot Delta kFS vs kB for different W and convert W into rms B1. 
% figure(3)
% for pp = 1:length(W)
%     plot(k_B, k_FSapp(pp,:),'--','LineWidth', 2)
%     hold on
% end
% 
% xlabel('k_{B} [s^{-1}]', 'FontSize', 18)
% ylabel('\Deltak_{FS} [s^{-1}]', 'FontSize', 18) 
% l = legend('0.792', '1.77', '2.16', '2.50', '3.07');
% title(l, 'B_{1,rms} [\muT]','FontSize', 18)
% l.FontSize = 18;
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)

figure(2)
cm = colormap(cool(length(B1)));
plot(B1rms*1e6, T1_Fapp,'LineWidth', 2, 'color',cm(1,:))
hold on
grid on
grid minor
plot(B1rms*1e6, T1_Sapp,'LineWidth', 2, 'color',cm(29,:))
xlabel('B_{1,rms} [\muT]', 'FontSize', 18)
ylabel('T_{1}^{app} [s]', 'FontSize', 18) 
l = legend('T_{1F}^{app}', 'T_{1S}^{app}');
xlim([B1rms(1)*1e6 B1rms(29)*1e6])
l.FontSize = 18;
get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)

% subplot(1,2,2)
% for pp = 1:length(W)
%     plot(k_B, T1_Sapp(pp,:),'LineWidth', 2, 'color',cm(pp,:))
%     hold on
% end
% xlabel('k_{B} [s^{-1}]', 'FontSize', 18)
% ylabel('T_{1S}^{app} [s]', 'FontSize', 18) 
% l = legend('0.792', '1.77', '2.16', '2.50', '3.07');
% title(l, 'B_{1,rms} [\muT]','FontSize', 18)
% l.FontSize = 18;
% get(gca, 'XTick'); get(gca, 'YTick'); set(gca, 'FontSize', 16)
