T1_F = 0.45; T1_S = 1.4; T2_F = 0.015; T2_S = 0.09; M0_F = 0.15; k_FS = 0:1:20; Delta = 0; PC1 = 0 + Delta; PC2 = pi + Delta;

cm = colormap(magma(5));
 
subplot(2,3,1); plot(k_FS, Solution_B(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1_S, 'g--'); get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14); grid on; grid minor;
subplot(2,3,2); plot(k_FS, Solution_B(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1_F, 'g--'); get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14); grid on; grid minor;
subplot(2,3,3); plot(k_FS, Solution_B(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2_S, 'g--'); get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14); grid on; grid minor;
subplot(2,3,4); plot(k_FS, Solution_B(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2_F, 'g--'); get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14); grid on; grid minor;
subplot(2,3,5); plot(k_FS, Solution_B(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('MWF', 'FontSize', 12); hline(M0_F, 'g--'); get(gca, 'XTick'); set(gca, 'FontSize', 14); get(gca, 'YTick'); set(gca, 'FontSize', 14); grid on; grid minor;
%subplot(2,3,6); plot(k_FS, Solution_NE_Bouhrara(:,:,6), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('k_{FS} (s^{-1})', 'FontSize', 12);
 
subplot(2,3,1); plot(k_FS, Solution_D(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_D(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,3); plot(k_FS, Solution_D(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,4); plot(k_FS, Solution_D(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_D(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
%subplot(2,3,6); plot(k_FS, Solution_NE_Deoni(:,:,6), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
 
subplot(2,3,1); plot(k_FS, Solution_W(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_W(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,3); plot(k_FS, Solution_W(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,4); plot(k_FS, Solution_W(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_W(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
%subplot(2,3,6); plot(k_FS, Solution_NE_Wood(:,:,6), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
 
subplot(2,3,1); plot(k_FS, Solution_Z(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,2); plot(k_FS, Solution_Z(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,3); plot(k_FS, Solution_Z(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,4); plot(k_FS, Solution_Z(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,5); plot(k_FS, Solution_Z(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
%subplot(2,3,6); plot(k_FS, Solution_NE_Zhang(:,:,6), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
ll = legend('B1','B2','B3','B4'); ll.FontSize = 16;
