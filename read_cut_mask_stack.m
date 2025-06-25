%% 读取、裁剪文件--------------------可调 folder_path（文件名），wavenumber（波数），row_range、col_range（裁剪范围）和threshold（阈值）
% 初始设置
clear; clc;
folder_path = '/Volumes/Public/Zhan Zhihao/20250624 SRS algae/843.26/04PA'; 
num_channels = 92;  % 通道数
wavenumbers = linspace(2000, 2300, num_channels); % 波数
row_range = 322:374;   % ROI: 行(Y)
col_range = 110:180;   % ROI: 列(X)
% 取一张图确认尺寸
sample_img = imread(fullfile(folder_path, '0.tif'));
roi_img = sample_img(row_range, col_range);
height = numel(row_range);
width  = numel(col_range);
stack = zeros(height, width, num_channels, 'like', sample_img);
stack = double(stack);

% 循环读取 + 裁剪
for ch = 0:num_channels-1
    filename = fullfile(folder_path, sprintf('%d.tif', ch));
    if isfile(filename)
        img = imread(filename);
        roi_img = img(row_range, col_range);
        stack(:,:,ch+1) = roi_img;
    else
        warning('缺失文件: %s', filename);
    end
end

% 像素总数
pixel_amount = width * height;

%% 阈值分割及后处理
% 选择参考通道图像（用于分割藻细胞区域）
ref_channel = 39;
ref_img = double(stack(:, :, ref_channel));  % 转换为 double 类型以便计算
% 初始化邻域平均图像
[height, width] = size(ref_img);
avg_img = zeros(height, width);
% 补零边界处理，避免越界
ref_padded = padarray(ref_img, [1, 1], 'replicate');
% 计算每个像素与上下左右共 5 像素的平均值
for i = 2:height+1
    for j = 2:width+1
        % 中心 + 上下左右
        neighborhood = [ ...
            ref_padded(i, j), ...
            ref_padded(i-1, j), ...
            ref_padded(i+1, j), ...
            ref_padded(i, j-1), ...
            ref_padded(i, j+1) ...
        ];
        avg_img(i-1, j-1) = mean(neighborhood);
    end
end
% 手动设定阈值（可调）
threshold = 32900;
% 平均强度高于阈值认为是藻细胞区域
cell_mask = avg_img > threshold;
% 后处理 平滑形状、填充孔洞、去小噪声 
se = strel('disk', 2);
cell_mask_clean = imopen(cell_mask, se);
cell_mask_clean = imfill(cell_mask_clean, 'holes');
cell_mask_clean = bwareaopen(cell_mask_clean, 50);
figure; imshow(cell_mask_clean, [],'InitialMagnification', 500);
title('Refined Cell Mask' );

% % 不进行后处理直接成图（对比）
% figure; imshow(cell_mask, [],'InitialMagnification', 500);
% title('Unrefined Cell Mask' );

%% 算背景平均谱、扣除背景平均谱
% 拉直成二维：[像素数 × 通道数]
stack_reshaped = reshape(stack, [], num_channels);  % [H*W, 92]
%把分割后的 mask 展平成列向量
mask_flat = cell_mask_clean(:);  % [N_pixels × 1]
% mask_flat = cell_mask (:);
% 提取背景
bg_spectra = stack_reshaped(~mask_flat, :);  % ~mask_flat 表示 mask=0 的部分
% 求背景平均谱
background_spectrum = mean(bg_spectra, 1);  % [1 × num_channels]
% % 可视化
% figure;
% plot(wavenumbers, background_spectrum, 'LineWidth', 1.8);
% xlabel('Raman Shift (cm^{-1})');
% ylabel('Intensity');
% title('Mean Spectrum of Background');
% grid on;
% 对 stack 的每个像素谱线扣去背景谱
corrected_stack = zeros(size(stack));
for c = 1:num_channels
    corrected_stack(:, :, c) = stack(:, :, c) - background_spectrum(c);    %这里之后用corrected_stack来表示细胞
end

%% 自动保存文件名
[~, folder_name] = fileparts(folder_path);
save_name = sprintf('%s_cut_corrected_stack.mat', folder_name);
save(fullfile(folder_path, save_name), 'stack', 'wavenumbers', 'cell_mask');
fprintf('已保存: %s\n', fullfile(folder_path, save_name));