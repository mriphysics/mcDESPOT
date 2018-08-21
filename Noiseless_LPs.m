cm = colormap(cool(4));
 
subplot(2,3,1); plot(k_FS, Solution_NE_Bouhrara(:,:,1), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1S} (s)', 'FontSize', 12); hline(T1_S, 'r')
subplot(2,3,2); plot(k_FS, Solution_NE_Bouhrara(:,:,2), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{1F} (s)', 'FontSize', 12); hline(T1_F, 'r')
subplot(2,3,3); plot(k_FS, Solution_NE_Bouhrara(:,:,3), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2S} (s)', 'FontSize', 12); hline(T2_S, 'r')
subplot(2,3,4); plot(k_FS, Solution_NE_Bouhrara(:,:,4), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('T_{2F} (s)', 'FontSize', 12); hline(T2_F, 'r')
subplot(2,3,5); plot(k_FS, Solution_NE_Bouhrara(:,:,5), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('M_{0F}', 'FontSize', 12); hline(M0_F, 'r')
%subplot(2,3,6); plot(k_FS, Solution_NE_Bouhrara(:,:,6), '--o', 'Color',  cm(1,:), 'LineWidth', 2); hold on;  xlabel('GT k_{FS} (s^{-1})', 'FontSize', 12); ylabel('k_{FS} (s^{-1})', 'FontSize', 12);
 
subplot(2,3,1); plot(k_FS, Solution_NE_Deoni(:,:,1), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_NE_Deoni(:,:,2), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,3); plot(k_FS, Solution_NE_Deoni(:,:,3), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on; 
subplot(2,3,4); plot(k_FS, Solution_NE_Deoni(:,:,4), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_NE_Deoni(:,:,5), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
%subplot(2,3,6); plot(k_FS, Solution_NE_Deoni(:,:,6), '--o', 'Color',  cm(2,:), 'LineWidth', 2); hold on;  
 
subplot(2,3,1); plot(k_FS, Solution_NE_Wood(:,:,1), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on; 
subplot(2,3,2); plot(k_FS, Solution_NE_Wood(:,:,2), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,3); plot(k_FS, Solution_NE_Wood(:,:,3), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,4); plot(k_FS, Solution_NE_Wood(:,:,4), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
subplot(2,3,5); plot(k_FS, Solution_NE_Wood(:,:,5), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
%subplot(2,3,6); plot(k_FS, Solution_NE_Wood(:,:,6), '--o', 'Color',  cm(3,:), 'LineWidth', 2); hold on;  
 
subplot(2,3,1); plot(k_FS, Solution_NE_Zhang(:,:,1), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,2); plot(k_FS, Solution_NE_Zhang(:,:,2), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,3); plot(k_FS, Solution_NE_Zhang(:,:,3), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,4); plot(k_FS, Solution_NE_Zhang(:,:,4), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
subplot(2,3,5); plot(k_FS, Solution_NE_Zhang(:,:,5), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
%subplot(2,3,6); plot(k_FS, Solution_NE_Zhang(:,:,6), '--o', 'Color',  cm(4,:), 'LineWidth', 2);
ll = legend('Set 1','Set 2','Set 3','Set 4'); ll.FontSize = 12;
