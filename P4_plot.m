% 清除工作区和图形窗口
clear; clc; close all;

% --- 数据定义 ---

% 信噪比 SNR (dB)
snr = -4 : 0.5 : 0.5;

% 方案 1: Alamouti (Uncoded) 数据
ber_alamouti = [1.21e-01, 1.09e-01, 9.55e-02, 8.51e-02, 7.60e-02, ...
                6.48e-02, 5.61e-02, 4.88e-02, 4.10e-02, 3.29e-02];

% 方案 2: SM + 5G LDPC (K=520, Rate=0.5, BG2) 数据
ber_ldpc_520 = [3.10e-01, 2.94e-01, 2.91e-01, 2.62e-01, 2.46e-01, ...
                1.88e-01, 1.17e-01, 4.51e-02, 6.75e-03, 5.70e-04];

% 方案 3: SM + 5G LDPC (K=1040, Rate=0.5, BG2) 数据
ber_ldpc_1040 = [3.13e-01, 2.94e-01, 2.75e-01, 2.66e-01, 2.45e-01, ...
                 2.12e-01, 1.07e-01, 2.48e-02, 6.46e-04, 8.17e-06];

% --- 绘图部分 ---

figure;
% 绘制 Alamouti 曲线 (黑色方块)
semilogy(snr, ber_alamouti, '-ks', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
hold on;

% 绘制 LDPC K=520 曲线 (蓝色圆圈)
semilogy(snr, ber_ldpc_520, '-bo', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'b');

% 绘制 LDPC K=1040 曲线 (红色菱形)
semilogy(snr, ber_ldpc_1040, '-rd', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'r');

% --- 图形美化 ---

grid on; % 打开网格
xlabel('SNR (dB)', 'FontSize', 5);
ylabel('Bit Error Rate (BER)', 'FontSize', 5);
title('BER Comparison: Alamouti vs. 5G LDPC', 'FontSize', 5);

% 设置图例
legend('Uncoded Alamouti', ...
       'SM + 5G LDPC (K=520, R=0.5)', ...
       'SM + 5G LDPC (K=1040, R=0.5)', ...
       'Location', 'SouthWest');

% 设置坐标轴范围 (根据数据微调以获得最佳视觉效果)
xlim([-4 0.5]);
ylim([1e-6 1]);

% 开启次要网格线，方便读数
grid minor;

set(gca, 'FontSize', 12);
hold off;