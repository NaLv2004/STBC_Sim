%% Massive MIMO vs Small MIMO: 信道硬化多角度仿真演示
% 包含三个部分：
% 1. 信道增益 PDF (展示方差收敛) -> Subplot 1
% 2. 容量 CDF (展示中断概率/可靠性) -> Subplot 2
% 3. 矩阵条件数 (展示线性检测可行性) -> Separate Figure
% -----------------------------------------------------------

clear; clc; close all;

%% ================== 通用绘图设置 ==================
lineWidth = 2;
fontSize = 12;
fontName = 'Arial';

fprintf('仿真开始，正在生成图表...\n');

% 创建一个宽屏 Figure 用于容纳两个子图
figure(1); 
set(gcf, 'Color', 'w', 'Position', [100, 100, 1200, 500]); 

%% ================== 子图 1: 归一化信道增益 PDF ==================
% 目标：展示随着天线数增加，信道能量如何集中到一个确定值
subplot(1, 2, 1); % <--- 左侧子图
hold on; grid on; box on;

fprintf('正在绘制子图 1: 信道增益 PDF...\n');

Nr_list = [4, 16, 64, 128]; % 接收天线数量列表
n_realizations = 10000;      % 蒙特卡洛次数
colors = {'#0072BD', '#D95319', '#EDB120', '#77AC30'}; % 专业的配色

for i = 1:length(Nr_list)
    Nr = Nr_list(i);
    % 生成信道 h ~ CN(0, 1)
    h = (randn(Nr, n_realizations) + 1j * randn(Nr, n_realizations)) / sqrt(2);
    
    % 硬化指标：归一化能量 ||h||^2 / Nr
    h_power_norm = sum(abs(h).^2, 1) / Nr;
    
    % 计算 PDF
    [f, x] = ksdensity(h_power_norm);
    
    % 绘图
    plot(x, f, 'LineWidth', lineWidth, 'Color', colors{i}, ...
        'DisplayName', ['N_r = ' num2str(Nr)]);
end

xlabel('Normalized Channel Gain $(\|\bf{h}\|\|^2 / N_r)$', 'FontSize', fontSize, 'FontName', fontName, 'Interpreter','latex');
ylabel('Probability Density (PDF)', 'FontSize', fontSize, 'FontName', fontName);
title('Channel Hardening: Variance Reduction', 'FontSize', 14, 'FontName', fontName);
legend('show', 'Location', 'NorthEast');
xlim([0.2 1.8]); 
ylim([0 8]); 

%% ================== 子图 2: 容量 CDF (可靠性对比) ==================
% 目标：展示 Massive MIMO 的容量曲线如何变成“阶跃函数”
subplot(1, 2, 2); % <--- 右侧子图
hold on; grid on; box on;

fprintf('正在绘制子图 2: 容量 CDF...\n');

Nt = 4;           % 用户数
SNR_dB = 10;      % 平均 SNR
SNR = 10^(SNR_dB/10);
num_trials = 5000;

configs = [4, 128]; % 对比 4x4 和 4x128
config_names = {'Small MIMO (4\times4)', 'Massive MIMO (4\times128)'};
line_styles = {'--', '-'};
col_cdf = {'b', 'r'};

for k = 1:length(configs)
    Nr = configs(k);
    capacities = zeros(1, num_trials);
    
    for i = 1:num_trials
        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2);
        
        % [关键点]: 归一化 H 以消除阵列增益，仅对比“分布形状”
        H_norm = H / sqrt(Nr/Nt); 
        
        % Shannon Capacity
        capacities(i) = real(log2(det(eye(Nt) + (SNR/Nt) * (H_norm' * H_norm))));
    end
    
    % 手动计算 CDF 以获得更好的绘图控制
    sorted_C = sort(capacities);
    y_cdf = (1:num_trials) / num_trials;
    
    plot(sorted_C, y_cdf, 'Color', col_cdf{k}, 'LineStyle', line_styles{k}, ...
        'LineWidth', lineWidth, 'DisplayName', config_names{k});
end

xlabel('Instantaneous Capacity (bps/Hz)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('CDF (Outage Probability)', 'FontSize', fontSize, 'FontName', fontName);
title('Reliability: Capacity Hardening', 'FontSize', 14, 'FontName', fontName);
legend('show', 'Location', 'NorthWest');

%% ================== 图 2: 矩阵条件数 (工程代价) ==================
% 目标：展示信道矩阵逐渐趋向正交，利于线性检测
% 这是一个单独的图，不包含在上面的 subplot 中
fprintf('正在绘制图 2: 矩阵条件数...\n');

K = 8;                  % 用户数
Nr_range = 8:8:128;     % 扫描天线数
num_trials_cond = 200;
cond_results = zeros(1, length(Nr_range));

for i = 1:length(Nr_range)
    Nr = Nr_range(i);
    temp_cond = zeros(1, num_trials_cond);
    for j = 1:num_trials_cond
        H = (randn(Nr, K) + 1j * randn(Nr, K)) / sqrt(2);
        temp_cond(j) = cond(H); % 计算条件数
    end
    cond_results(i) = mean(temp_cond); % 取平均
end

figure(2); set(gcf, 'Color', 'w', 'Position', [400, 100, 600, 400]);
semilogy(Nr_range ./ K, cond_results, '-o', 'Color', 'k', ...
    'MarkerFaceColor', 'r', 'LineWidth', lineWidth);
grid on; box on;

xlabel('Antenna Ratio (N_r / K)', 'FontSize', fontSize, 'FontName', fontName);
ylabel('Average Condition Number (Log Scale)', 'FontSize', fontSize, 'FontName', fontName);
title('Asymptotic Orthogonality', 'FontSize', 14, 'FontName', fontName);
yline(1, '--b', 'Perfect Orthogonality', 'LineWidth', 1.5);

fprintf('仿真完成！请查看生成的图表。\n');