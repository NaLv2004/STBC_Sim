% 清除工作区
clear; clc; close all;

% --- 数据定义 ---

% 1. 定义两组不同的 SNR 向量
snr_ldpc = -4 : 0.5 : 0.5;   % LDPC 数据的 SNR 轴
snr_polar = -4 : 1 : 3;      % Polar 和 Alamouti 数据的 SNR 轴

% 2. Alamouti 数据 (使用范围更广的那组数据)
ber_alamouti = [1.21e-01, 9.72e-02, 7.41e-02, 5.64e-02, 4.12e-02, ...
                2.76e-02, 1.79e-02, 1.14e-02];

% 3. 5G LDPC 数据
ber_ldpc_520  = [3.10e-01, 2.94e-01, 2.91e-01, 2.62e-01, 2.46e-01, ...
                 1.88e-01, 1.17e-01, 4.51e-02, 6.75e-03, 5.70e-04];
ber_ldpc_1040 = [3.13e-01, 2.94e-01, 2.75e-01, 2.66e-01, 2.45e-01, ...
                 2.12e-01, 1.07e-01, 2.48e-02, 6.46e-04, 8.17e-06];

% 4. 5G Polar 数据
ber_polar_64  = [4.40e-01, 3.89e-01, 2.49e-01, 1.28e-01, 3.25e-02, ...
                 6.97e-03, 1.47e-03, 1.65e-04];
ber_polar_128 = [4.81e-01, 4.16e-01, 3.38e-01, 1.31e-01, 2.67e-02, ...
                 3.32e-03, 4.24e-04, 4.61e-05];

ber_STTC = [0.374, 0.266, 0.149, 0.0596, 0.0169, 0.00347, 0.000750, 0.000204];

% --- 绘图部分 ---

figure;

% 1. 绘制 Alamouti (黑色方块) - 基准
semilogy(snr_polar, ber_alamouti, '-ks', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
hold on;

semilogy(snr_polar, ber_STTC, '-ko', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'k');

% 2. 绘制 LDPC (蓝色)
% K=520 (圆圈)
semilogy(snr_ldpc, ber_ldpc_520, '-bo', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'b');
% K=1040 (菱形)
semilogy(snr_ldpc, ber_ldpc_1040, '-bd', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'b');

% 3. 绘制 Polar (红色)
% K=64 (圆圈)
semilogy(snr_polar, ber_polar_64, '-ro', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'r');
% K=128 (菱形)
semilogy(snr_polar, ber_polar_128, '-rd', ...
    'LineWidth', 1.3, 'MarkerSize', 3, 'MarkerFaceColor', 'r');

% --- 样式设置 (严格遵守字号5、不加粗) --

grid on;
grid minor;

% 坐标轴标签
xlabel('SNR (dB)', 'FontSize', 9, 'FontWeight', 'normal');
ylabel('Bit Error Rate (BER)', 'FontSize', 9, 'FontWeight', 'normal');

% 标题
title('Performance Comparison: Alamouti vs. STTC vs. Polar', ...
      'FontSize', 10, 'FontWeight', 'normal');

% 图例
lgd = legend(...
    'Uncoded Alamouti', ...
    'Uncoded STTC',...
    'LDPC (K=520, R=0.5)', ...
    'LDPC (K=1040, R=0.5)', ...
    'Polar (K=64, E=128)', ...
    'Polar (K=128, E=256)', ...
    'Location', 'SouthWest');
set(lgd, 'FontSize', 5);

% 坐标轴刻度字体
%set(gca, 'FontSize', 5, 'FontWeight', 'normal');

% 调整范围
xlim([-4 3]);
ylim([1e-6 1]); % 调整下限以容纳LDPC最低的误码率

hold off;