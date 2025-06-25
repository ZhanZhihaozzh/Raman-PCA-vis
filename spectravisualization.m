clear; clc;
%% 读取
% 设置文件名
folder_path = '/Volumes/Public/Zhan Zhihao/20250624 SRS algae/843.26/04PA';
% 自动提取保存名
[~, folder_name] = fileparts(folder_path);
save_name = sprintf('%s_cut_corrected_stack.mat', folder_name);
data = load(fullfile(folder_path, save_name));
stack = data.stack;
wavenumbers = data.wavenumbers;
cell_mask = data.cell_mask;
% 
% % 显示其中一个通道（如第30个）
% figure;
% imshow(stack(:,:,38), []);
% title('1527cm^-^1 Channel 38');
% 
% % 提取并显示单一背景像素的光谱
% x = 100; y = 100;
% spectrum = squeeze(stack(y, x, :));
% figure;
% plot(wavenumbers,spectrum);
% xlabel('Raman Shift (cm^{-1})'); ylabel('SRS Intensity');
% title(sprintf('Spectrum at pixel (%d, %d)', x, y));
% grid on;
% 
% % 提取并显示单一微藻像素的光谱
% x = 100; y = 100;
% spectrum = squeeze(stack(y, x, :));
% figure;
% plot(wavenumbers,spectrum);
% xlabel('Raman Shift (cm^{-1})'); ylabel('SRS Intensity');
% title(sprintf('Spectrum at pixel (%d, %d)', x, y));
% grid on;

%% 可视化背景平均谱
% 求背景区域平均谱线
bg_mask = ~cell_mask;  % 背景mask
% 拉直成二维：[像素数 × 通道数]
[height, width, num_channels] = size(stack);
stack_reshaped = reshape(stack, [], num_channels);  % [H*W, 92]
mask_flat = reshape(bg_mask, [], 1);  % [H*W, 1]
% 找到背景像素
bg_spectra = stack_reshaped(mask_flat, :);  % [N_bg, 92]
background_spectrum = mean(bg_spectra, 1);  % [1 × 92]
% 可视化背景谱
figure;
plot(wavenumbers,background_spectrum, 'LineWidth', 1.8);
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity');
title('Mean Spectrum of Background');
grid on;
% 
% % 绘制某微藻像素点背景扣除后的谱线
% x = 148; y = 265;  
% corrected_spectrum = squeeze(corrected_stack(y, x, :));
% 
% figure;
% plot(corrected_spectrum, 'LineWidth', 1.5);
% xlabel('Channel');
% ylabel('Background-corrected Intensity');
% title(sprintf('Corrected Spectrum at (%d, %d)', x, y));

%% 可视化微藻所有像素去背景平均之后的平均谱线
% 强制转换为 double 类型（防止整数类型运算错误）
cell_mask_flat = reshape(cell_mask, [], 1)
stack_reshaped = reshape(stack, [], 92)
cell_spectra = stack_reshaped(cell_mask_flat, :)
cell_spectra = double(cell_spectra);
background_spectrum = double(background_spectrum);
% 可视化细胞谱
figure;
cell_mean = mean(cell_spectra,1);
plot(wavenumbers,cell_mean, 'LineWidth', 1.8);
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity');
title('Mean Spectrum of Cell');
grid on;
% 背景扣除
cell_corrected = cell_spectra - background_spectrum;  % [N_cell × 92]、
%平均
cell_corrected_mean = mean(cell_corrected, 1);
%绘图
figure;
plot(wavenumbers,cell_corrected_mean, 'LineWidth', 1.8);
xlabel('Raman Shift (cm^{-1})');
ylabel('Intensity (Background-subtracted)');
title('Mean Spectrum of Microalgae Cells (Background Subtracted)');
grid on;

% %% =============== 2. 先对全图做 PCA ================
% [coeff_all, score_all, latent_all, ~, explained_all] = pca(corrected_stack);
% H = 500;
% W = 500;
% % 方差贡献率
% figure;
% plot(cumsum(explained_all), '-o', 'LineWidth', 1.5);
% xlabel('Number of Principal Components');
% ylabel('Cumulative Variance Explained (%)');
% title('PCA - Whole Image: Cumulative Explained Variance');
% grid on;
% 
% % PC1 score 热图
% PC1_img_all = reshape(score_all(:,1), H, W);
% figure;
% imagesc(PC1_img_all);
% axis image off;
% colorbar;
% title('PCA - Whole Image: PC1 Score Map');
% 
% % PC2 score 热图
% PC2_img_all = reshape(score_all(:,2), H, W);
% figure;
% imagesc(PC2_img_all);
% axis image off;
% colorbar;
% title('PCA - Whole Image: PC2 Score Map');
% 
% % PC1 vs PC2 scatter plot
% figure;
% scatter(score_all(:,1), score_all(:,2), 5, '.');
% xlabel('PC1'); ylabel('PC2');
% title('PCA - Whole Image: PC1 vs PC2 Scatter');
% grid on;
% 
% %% =============== 4. 对细胞像素做 PCA ===============
% % % 以通道30 + 5像素邻域均值为例（快速版）
% % ref_channel = 30;
% % ref_img = stack(:,:,ref_channel);
% % 
% % kernel = [0 1 0; 1 1 1; 0 1 0]/5;
% % avg_img = conv2(ref_img, kernel, 'same');
% % 
% % threshold = 100;  % 请根据实际信号调节
% % cell_mask = avg_img > threshold;
% % 
% % % 展平 mask
% % cell_mask_flat = reshape(cell_mask, [], 1);
% % 
% % % 可视化分割效果
% % figure;
% % imshow(cell_mask);
% % title('Cell Mask');
% 
% % 做 PCA
% [coeff_cell, score_cell, latent_cell, ~, explained_cell] = pca(cell_spectra);
% 
% % 方差贡献率
% figure;
% plot(cumsum(explained_cell), '-o', 'LineWidth', 1.5);
% xlabel('Number of Principal Components');
% ylabel('Cumulative Variance Explained (%)');
% title('PCA - Cells Only: Cumulative Explained Variance');
% grid on;
% 
% % PC1 score 热图（只在细胞区域上显示）
% PC1_cell_img = nan(N, 1);
% PC1_cell_img(cell_mask_flat) = score_cell(:,1);
% PC1_cell_img = reshape(PC1_cell_img, H, W);
% 
% figure;
% imagesc(PC1_cell_img);
% axis image off;
% colorbar;
% title('PCA - Cells Only: PC1 Score Map');
% 
% % PC2 score 热图（只在细胞区域上显示）
% PC2_cell_img = nan(N, 1);
% PC2_cell_img(cell_mask_flat) = score_cell(:,2);
% PC2_cell_img = reshape(PC2_cell_img, H, W);
% 
% % PC1 vs PC2 vs PC3 (3D scatter) 
% figure;
% scatter3(score_cell(:,1), score_cell(:,2), score_cell(:,3), 8, '.');
% xlabel('PC1');
% ylabel('PC2');
% zlabel('PC3');
% title('PCA - Cells Only: PC1 vs PC2 vs PC3 (3D Scatter)');
% grid on;
% axis tight; 
% 
% figure;
% imagesc(PC2_cell_img);
% axis image off;
% colorbar;
% title('PCA - Cells Only: PC2 Score Map');
% 
% % PC1 vs PC2 scatter plot (cells only)
% figure;
% scatter(score_cell(:,1), score_cell(:,2), 8, '.');
% xlabel('PC1'); ylabel('PC2');
% title('PCA - Cells Only: PC1 vs PC2 Scatter');
% grid on;
