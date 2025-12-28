% =========================================================================
% Refined MIMO Simulation
% Scenario: 2x2 MIMO
% Part 1: Fair comparison at Fixed Rate (2 bits/s/Hz)
% Part 2: Asymptotic DMT Tradeoff (Alamouti vs Optimal 2x2)
% =========================================================================

clear; clc; close all;

%% =========================================================================
%% PART 1: SER vs SNR (Fair Comparison: Same Spectral Efficiency)
%% =========================================================================
fprintf('Running Part 1: Fair SER comparison (Rate = 2 bits/s/Hz)...\n');

SNR_dB_vec = 0:2:10;
SNR_lin_vec = 10.^(SNR_dB_vec./10);
num_trials = 2000; 

qpsk_syms = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2); 
SM_Codebook = [1, 1, -1, -1; 1, -1, 1, -1]; % BPSK combinations

ser_sm = zeros(size(SNR_dB_vec));
ser_ala = zeros(size(SNR_dB_vec));

for s_idx = 1:length(SNR_dB_vec)
    snr_lin = SNR_lin_vec(s_idx);
    scale = sqrt(snr_lin / 2);
    
    err_sm = 0; err_ala = 0;
    
    for t = 1:num_trials
        H = (randn(2,2) + 1i*randn(2,2)) / sqrt(2);
        n = (randn(2,1) + 1i*randn(2,1)) / sqrt(2);
        
        % -- SM (BPSK) --
        tx_idx = randi(4);
        y_sm = scale * H * SM_Codebook(:, tx_idx) + n;
        [~, est_idx] = min(sum(abs(y_sm - scale * H * SM_Codebook).^2, 1));
        if est_idx ~= tx_idx, err_sm = err_sm + 1; end
        
        % -- Alamouti (QPSK) --
        idx1 = randi(4); s1 = qpsk_syms(idx1);
        h_pow = norm(H, 'fro')^2;
        n_eff = (randn + 1i*randn)/sqrt(2);
        y_ala = sqrt(scale^2 * h_pow) * s1 + n_eff;
        [~, det_i] = min(abs(y_ala - sqrt(scale^2 * h_pow) * qpsk_syms));
        if det_i ~= idx1, err_ala = err_ala + 1; end
    end
    ser_sm(s_idx) = err_sm / num_trials;
    ser_ala(s_idx) = err_ala / num_trials;
end

%% =========================================================================
%% PART 2: Diversity-Multiplexing Tradeoff (Semi-Analytic)
%% =========================================================================
fprintf('\nRunning Part 2: Asymptotic DMT Tradeoff...\n');

SNR_dB_DMT = [200, 210, 220]; 
SNR_lin_DMT = 10.^(SNR_dB_DMT./10);

% Alamouti simulation range (r from 0 to 1)
r_vals = 0:0.1:1; 
d_simulated = zeros(size(r_vals));

for r_idx = 1:length(r_vals)
    r = r_vals(r_idx);
    
    % Fix for r=0 calculation (avoid log(0) or 0 slope)
    r_calc = max(r, 1e-10); 
    
    p_out = zeros(size(SNR_dB_DMT));
    for s_i = 1:length(SNR_dB_DMT)
        snr = SNR_lin_DMT(s_i);
        % Threshold for Alamouti: X < 2 * (SNR^r - 1) / SNR
        threshold_X = 2 * (snr^r_calc - 1) / snr;
        % CDF of ||H||_F^2 for 2x2 (Gamma(4,1))
        p_out(s_i) = gamcdf(threshold_X, 4, 1); 
        if p_out(s_i) == 0, p_out(s_i) = 1e-100; end
    end
    
    coeffs = polyfit(log10(SNR_lin_DMT), log10(p_out), 1);
    d_simulated(r_idx) = -coeffs(1);
end

% --- Theoretical Curves ---
% 1. Alamouti Theory
d_ala_theory = 4 * (1 - r_vals);

% 2. Optimal 2x2 MIMO Theory (Zheng & Tse)
% Points: (0,4) -> (1,1) -> (2,0)
r_optimal = [0, 1, 2];
d_optimal = [4, 1, 0];

%% =========================================================================
%% PLOTTING
%% =========================================================================

figure('Position', [100, 100, 1200, 500]);

% % Plot 1: SER
% subplot(1, 2, 1);
% semilogy(SNR_dB_vec, ser_sm, 'b-o', 'LineWidth', 2); hold on;
% semilogy(SNR_dB_vec, ser_ala, 'r-s', 'LineWidth', 2);
% grid on;
% xlabel('SNR (dB)');
% ylabel('Symbol Error Rate (SER)');
% title({'Fair Comparison (Rate = 2 bits/s/Hz)', 'Alamouti(QPSK) vs SM(BPSK)'});
% legend('2x2 SM (BPSK)', '2x2 Alamouti (QPSK)', 'Location', 'SouthWest');

% Plot 2: DMT
subplot(1, 1, 1);

% A. Optimal Curve (Green)
plot(r_optimal, d_optimal, 'g-', 'LineWidth', 3); hold on;

% B. Alamouti Theory (Black Dashed)
plot(r_vals, d_ala_theory, 'k--', 'LineWidth', 2); 

% C. Alamouti Simulated (Blue Dots)
plot(r_vals, d_simulated, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');

grid on;
xlabel('Multiplexing Gain r');
ylabel('Diversity Gain d');
title({'Diversity-Multiplexing Tradeoff', '2x2 MIMO System'});

% Adjust axis to show full 2x2 capability (r up to 2)
axis([0 2 0 4.5]); 
xticks(0:0.2:2);

legend('Optimal 2x2 Tradeoff (Zheng & Tse)', ...
       'Alamouti Theory d(r)=4(1-r)', ...
       'Alamouti Simulated', ...
       'Location', 'NorthEast');