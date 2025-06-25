clear; clc;
%% 读取
% 设置文件名
folder_path = '/Volumes/Public/Zhan Zhihao/20250624 SRS algae/843.26/03PA';
% 自动提取保存名
[~, folder_name] = fileparts(folder_path);
save_name = sprintf('%s_cut_corrected_stack.mat', folder_name);
data = load(fullfile(folder_path, save_name));
stack = data.stack;
wavenumbers = data.wavenumbers;
cell_mask = data.cell_mask;

%% PCA--All field
%展平 
[H, W, C] = size(stack);
N = H * W;
stack_reshaped = reshape(stack, [], C);
% PCA 
[coeff, score, latent, ~, explained] = pca(stack_reshaped);
% 贡献率
figure;
subplot(1,2,1);
bar(explained);
xlabel('PC Index');
ylabel('Variance Explained (%)');
title('Individual PC Contribution');

subplot(1,2,2);
plot(cumsum(explained), '-o');
xlabel('Number of PCs');
ylabel('Cumulative Variance Explained (%)');
title('Cumulative Contribution');
grid on;
%PC 图像
for i = 1:3
    PC_map = reshape(score(:,i), H, W);
    figure;
    imagesc(PC_map);
    axis image off; colorbar;
    title(sprintf('PC%d Image', i));
end
% PC scatter plot 
figure;
scatter(score(:,1), score(:,2), 5, '.');
xlabel('PC1'); ylabel('PC2');
title('PC1 vs PC2 Scatter');

%% PCA--cell
% 展平堆栈与 mask
[H, W, C] = size(stack);
N = H * W;
stack_reshaped = reshape(stack, N, C);
mask_flat = cell_mask(:);
% 提取细胞像素
cell_spectra = stack_reshaped(mask_flat, :);  % [N_cells × C]

% PCA 
[coeff, score, latent, ~, explained] = pca(cell_spectra);

% 贡献率
figure;
subplot(1,2,1);
bar(explained);
xlabel('PC Index');
ylabel('Variance Explained (%)');
title('Individual PC Contribution');
subplot(1,2,2);
plot(cumsum(explained), '-o');
xlabel('Number of PCs');
ylabel('Cumulative Variance Explained (%)');
title('Cumulative Contribution');
grid on;

%PC图像
for i = 1:3
figure;
PC_map_flat = nan(H*W, 1);  % 改成 NaN
PC_map_flat(mask_flat) = score(:, i);
PC_map = reshape(PC_map_flat, H, W);
imagesc(PC_map);
axis image off; colorbar;
title(sprintf('PC%d Image (NaN background)', i));
end

% scatter 3D
pc1 = score(:,1);
pc2 = score(:,2);
pc3 = score(:,3);
range1 = max(pc1) - min(pc1);
range2 = max(pc2) - min(pc2);
range3 = max(pc3) - min(pc3);
max_range = max([range1, range2, range3]);
mid1 = (max(pc1) + min(pc1)) / 2;
mid2 = (max(pc2) + min(pc2)) / 2;
mid3 = (max(pc3) + min(pc3)) / 2;
figure;
scatter3(pc1, pc2, pc3, 10, 'filled');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('3D Scatter: PC1 vs PC2 vs PC3 (Cells Only)');
grid on;
xlim([mid1 - max_range/2, mid1 + max_range/2]);
ylim([mid2 - max_range/2, mid2 + max_range/2]);
zlim([mid3 - max_range/2, mid3 + max_range/2]);
axis equal;
axis vis3d;
