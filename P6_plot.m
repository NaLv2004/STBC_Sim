function plot_mimo_subplots()
    % 清除环境
    clear; clc; close all;

    % =====================================================================
    % 1. 数据录入
    % =====================================================================
    snr_axis = -8:1:4;
    
    % --- T_coh = 2 (Fast Fading) ---
    data.T2.Ala_K64  = [4.47e-01, 4.30e-01, 3.31e-01, 1.66e-01, 5.50e-02, 1.20e-02, 1.87e-03, 2.09e-04, 5.00e-05, 1.56e-06, 0, 0, 0];
    data.T2.SM_K64   = [3.96e-01, 3.15e-01, 2.21e-01, 9.43e-02, 1.53e-02, 2.66e-03, 2.12e-04, 1.87e-05, 0, 0, 0, 0, 0];
    data.T2.Ala_K128 = [4.48e-01, 4.63e-01, 3.65e-01, 2.04e-01, 6.67e-02, 4.92e-03, 5.98e-04, 5.31e-05, 2.50e-05, 0, 0, 0, 0];
    data.T2.SM_K128  = [4.77e-01, 4.02e-01, 2.87e-01, 8.84e-02, 7.16e-03, 4.55e-04, 2.34e-05, 0, 0, 0, 0, 0, 0];

    % --- T_coh = 16 (Block Fading) ---
    data.T16.Ala_K64  = [4.66e-01, 3.82e-01, 2.97e-01, 1.76e-01, 6.45e-02, 1.87e-02, 5.55e-03, 1.02e-03, 1.09e-04, 3.13e-05, 0, 0, 0];
    data.T16.SM_K64   = [3.95e-01, 3.33e-01, 1.98e-01, 7.16e-02, 3.45e-02, 5.98e-03, 9.41e-04, 1.52e-04, 4.06e-05, 0, 0, 0, 0];
    data.T16.Ala_K128 = [4.83e-01, 4.14e-01, 3.67e-01, 1.81e-01, 5.46e-02, 1.10e-02, 1.42e-03, 1.38e-04, 4.22e-05, 0, 0, 0, 0];
    data.T16.SM_K128  = [4.46e-01, 3.76e-01, 2.47e-01, 9.19e-02, 1.67e-02, 1.23e-03, 5.47e-05, 0, 0, 0, 0, 0, 0];

    % --- T_coh = 64 (Slow Fading) ---
    data.T64.Ala_K64  = [3.92e-01, 3.78e-01, 3.01e-01, 1.65e-01, 9.04e-02, 5.98e-02, 2.37e-02, 9.16e-03, 3.96e-03, 9.53e-04, 4.31e-04, 2.81e-05, 0];
    data.T64.SM_K64   = [4.09e-01, 2.97e-01, 2.12e-01, 1.24e-01, 5.83e-02, 2.56e-02, 1.40e-02, 4.31e-03, 1.20e-03, 2.95e-04, 1.11e-04, 6.25e-06, 0];
    data.T64.Ala_K128 = [4.37e-01, 3.94e-01, 3.36e-01, 2.29e-01, 9.01e-02, 4.42e-02, 7.78e-03, 1.53e-03, 3.38e-04, 2.03e-05, 0, 0, 0];
    data.T64.SM_K128  = [4.48e-01, 4.02e-01, 2.18e-01, 1.44e-01, 6.46e-02, 1.49e-02, 3.31e-03, 3.93e-04, 1.64e-05, 0, 0, 0, 0];

    % =====================================================================
    % 2. 绘图设置
    % =====================================================================
    figure('Position', [50, 100, 1400, 500], 'Color', 'w');
    
    % 定义场景列表方便循环
    scenarios = {data.T2, data.T16, data.T64};
    titles = {'Fast Fading (T_{coh}=2)', 'Block Fading (T_{coh}=16)', 'Slow Fading (T_{coh}=64)'};
    
    for i = 1:3
        subplot(1, 3, i);
        hold on;
        
        curr_data = scenarios{i};
        
        % 数据清洗 (0 -> NaN)
        y1 = clean(curr_data.Ala_K64);
        y2 = clean(curr_data.Ala_K128);
        y3 = clean(curr_data.SM_K64);
        y4 = clean(curr_data.SM_K128);
        
        % --- 绘制曲线 ---
        
        % 1. Alamouti K=64 (R=0.5 -> N=128)
        % 样式：黑色虚线 + 圆圈
        p1 = plot(snr_axis, y1, '--ko', 'LineWidth', 1.5, 'MarkerSize', 6);
        
        % 2. Alamouti K=128 (R=0.5 -> N=256)
        % 样式：黑色实线 + 方块
        p2 = plot(snr_axis, y2, '-ks', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        
        % 3. SM K=64 (R=0.25 -> N=256)
        % 样式：蓝色虚线 + 圆圈
        p3 = plot(snr_axis, y3, '--bo', 'LineWidth', 1.5, 'MarkerSize', 6);
        
        % 4. SM K=128 (R=0.25 -> N=512)
        % 样式：蓝色实线 + 方块
        p4 = plot(snr_axis, y4, '-bs', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        
        % --- 装饰 ---
        grid on; grid minor;
        xlabel('SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
        if i == 1
            ylabel('Bit Error Rate (BER)', 'FontSize', 11, 'FontWeight', 'bold');
        end
        title(titles{i}, 'FontSize', 12);
        
        % 强制 Log 坐标
        set(gca, 'YScale', 'log');
        ylim([1e-6, 1]);
        xlim([-8, 4]);
        
        % 图例 (Legend)
        legend([p1, p2, p3, p4], ...
               'Alamouti (R=0.5, N=128, K=64)', ...
               'Alamouti (R=0.5, N=256, K=128)', ...
               'SM (R=0.25, N=256, K=64)', ...
               'SM (R=0.25, N=512, K=128)', ...
               'Location', 'SouthWest', 'FontSize', 9);
    end
    
    % 添加总标题
    sgtitle('MIMO Polar Codes Performance: Alamouti vs. Spatial Multiplexing (SM)', 'FontSize', 14, 'FontWeight', 'bold');
end

function y = clean(x)
    y = x;
    y(y==0) = NaN;
end