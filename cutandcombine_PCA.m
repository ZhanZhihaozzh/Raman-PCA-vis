%% 加载两个 corrected_stack 文件 
% 设置文件名1
folder_path = '/Volumes/Public/Zhan Zhihao/20250607 SRS algae/790/49PA';
% 自动提取保存名1
[~, folder_name] = fileparts(folder_path);
save_name = sprintf('%s_cut_corrected_stack.mat', folder_name);
stack1 = load(fullfile(folder_path, save_name));
% 设置文件名2
folder_path = '/Volumes/Public/Zhan Zhihao/20250624 SRS algae/790/52PA';
% 自动提取保存名2
[~, folder_name] = fileparts(folder_path);
save_name = sprintf('%s_cut_corrected_stack.mat', folder_name);
stack2 = load(fullfile(folder_path, save_name)); 
 
% ====== 从两个stack结构中提取字段 ======
data1 = stack1.stack;
mask1 = stack1.cell_mask;
data2 = stack2.stack;
mask2 = stack2.cell_mask;
wavenumbers = stack1.wavenumbers;  % 公用的波数轴

% ====== 统计两个文件的细胞像素数 ======
num_cells1 = nnz(mask1);
num_cells2 = nnz(mask2);

% ====== 裁剪像素较多的文件，使两者细胞像素数相同 ======
if num_cells1 > num_cells2
    [H1, W1, ~] = size(data1);
    diff = num_cells1 - num_cells2;
    if H1 >= W1
        % 沿垂直方向裁剪下半边
        data1 = data1(1:floor(H1/2), :, :);
        mask1 = mask1(1:floor(H1/2), :);
    else
        % 沿水平方向裁剪右半边
        data1 = data1(:, 1:floor(W1/2), :);
        mask1 = mask1(:, 1:floor(W1/2));
    end
elseif num_cells2 > num_cells1
    [H2, W2, ~] = size(data2);
    diff = num_cells2 - num_cells1;
    if H2 >= W2
        % 沿垂直方向裁剪下半边
        data2 = data2(1:floor(H2/2), :, :);
        mask2 = mask2(1:floor(H2/2), :);
    else
        % 沿水平方向裁剪右半边
        data2 = data2(:, 1:floor(W2/2), :);
        mask2 = mask2(:, 1:floor(W2/2));
    end
end

% ====== 转换为double并展开 ======
data1 = double(data1);
data2 = double(data2);
[H1, W1, C] = size(data1);
[H2, W2, ~] = size(data2);

% ====== 展平图像并提取细胞区域像素 ======
flat1 = reshape(data1, [], C);
flat2 = reshape(data2, [], C);
mask1_flat = reshape(mask1, [], 1);
mask2_flat = reshape(mask2, [], 1);

cells1 = flat1(mask1_flat, :);  % [N1 x C]
cells2 = flat2(mask2_flat, :);  % [N2 x C]

% ====== 合并两个数据集并执行PCA ======
all_cells = [cells1; cells2];
[coeff, score, latent] = pca(all_cells);

% ====== 绘制方差贡献率图 ====== 
figure;
subplot(1,2,1);
explained = latent / sum(latent) * 100;
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

% ====== 绘制3D散点图（用颜色区分两个组） ======
N1 = size(cells1, 1);
N2 = size(cells2, 1);
score1 = score(1:N1, 1:3);
score2 = score(N1+1:end, 1:3);

figure;
scatter3(score1(:,1), score1(:,2), score1(:,3), 10, [0 0 0.8], 'filled'); hold on;
scatter3(score2(:,1), score2(:,2), score2(:,3), 10, [0.8 0 0], 'filled');
legend('Control Group','Isotope Group');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('3D PCA Scatter Plot (Cells Only)');
grid on;
axis equal;

% ====== PCA载荷图（反映主成分与原始波数的关系） ======
figure;
plot(wavenumbers, coeff(:,1:3), 'LineWidth', 1.5);
xlabel('Raman Shift (cm^{-1})'); ylabel('Loading Value');
title('PCA Loadings');
legend('PC1','PC2','PC3');
grid on;

% ====== PC图像：展示每个主成分在图像空间的分布（对齐后） ======
H = max(H1, H2);
W = W1 + W2;  % 左右拼接图像宽度

% 为每个PC绘制空间图像
for pc_idx = 1:3
    % 初始化总图像为NaN
    PC_map_full = nan(H * W, 1);

    % 准备数据1中的主成分图像
    score_map1 = nan(H1 * W1, 1);
    score_map1(mask1_flat) = score1(:, pc_idx);
    img1 = reshape(score_map1, H1, W1);
    pad1 = nan(H, W1);
    start_row1 = floor((H - H1)/2) + 1;
    pad1(start_row1:start_row1+H1-1, :) = img1;

    % 准备数据2中的主成分图像
    score_map2 = nan(H2 * W2, 1);
    score_map2(mask2_flat) = score2(:, pc_idx);
    img2 = reshape(score_map2, H2, W2);
    pad2 = nan(H, W2);
    start_row2 = floor((H - H2)/2) + 1;
    pad2(start_row2:start_row2+H2-1, :) = img2;

    % 拼接成左右结构图像
    PC_map = [pad1, pad2];

    % 绘图
    figure;
    imagesc(PC_map);
    axis image off; colorbar;
    title(sprintf('PC%d Image (Aligned, Cells Only)', pc_idx));

    % 添加组标签
    text(W1/2, H + 10, 'Control Group', 'HorizontalAlignment', 'center');
    text(W1 + W2/2, H + 10, 'Isotope Group', 'HorizontalAlignment', 'center');
end
