cm = colormap(cool(4));

subplot(2,3,1); histogram(Solution(:,1), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,2), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(T1_S, 'k--'); vline(0.2, 'r--'); vline(0.7, 'r--')
xlabel('T_{1S} (s)','FontSize',12); ylabel('Count','FontSize',12);

subplot(2,3,2); histogram(Solution(:,2), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,2), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(T1_F, 'k--'); vline(0.2, 'r--'); vline(0.7, 'r--')
xlabel('T_{1F} (s)','FontSize',12); ylabel('Count','FontSize',12);
 
subplot(2,3,3); histogram(Solution(:,3), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,3), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(T2_S, 'k--'); vline(0.06, 'r--'); vline(0.16, 'r--')
xlabel('T_{2S} (s)','FontSize',12); ylabel('Count','FontSize',12);

subplot(2,3,4); histogram(Solution(:,4), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,4), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(T2_F, 'k--'); vline(0.002, 'r--'); vline(0.04, 'r--')
xlabel('T_{2F} (s)','FontSize',12); ylabel('Count','FontSize',12);

subplot(2,3,5); histogram(Solution(:,5), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,5), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(M0_F, 'k--'); vline(0, 'r--'); vline(0.5, 'r--')
xlabel('M_{0F}','FontSize',12); ylabel('Count','FontSize',12);
 
subplot(2,3,6); histogram(Solution(:,6), 'EdgeColor', 'k', 'FaceColor', 'c', 'FaceAlpha', 0.8); hold on
histogram(Solution_Flat(:,6), 'EdgeColor', 'k', 'FaceColor', 'm', 'FaceAlpha', 0.8); 
vline(k_FS, 'k--'); vline(0.5, 'r--'); vline(20, 'r--')
xlabel('k_{FS} (s^{-1})','FontSize',12); ylabel('Count','FontSize',12);
 
ll = legend('SRC on Informative CF','SRC on Non-Informative CF'); ll.FontSize = 12;
