% =========================================================================
% DMT Visualization via SER Slope
% Scenario: 2x2 MIMO with Alamouti Scheme
% Methodology: Constellation Contraction to simulate Multiplexing Gain r
% Note: Removed dependency on Statistics Toolbox
% =========================================================================

clear; clc; close all;

fprintf('Running DMT Slope Simulation...\n');
fprintf('Simulating QPSK with distance scaling: d ~ SNR^(-r/2)\n');

%% Parameters
SNR_dB = 0:2:20; % High SNR range to see the slope clearly
SNR_lin = 10.^(SNR_dB./10);

% Define different Multiplexing Gains to test
r_values = [0, 0.25, 0.5, 0.75, 1.0]; 

num_trials = 500000; % Monte Carlo trials

% QPSK Constellation (Unit Energy)
qpsk_syms = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2);

%% Simulation Loop
% Store results for plotting
results_sim = zeros(length(r_values), length(SNR_dB));
results_theory = zeros(length(r_values), length(SNR_dB));

for r_idx = 1:length(r_values)
    r = r_values(r_idx);
    fprintf('Simulating for Multiplexing Gain r = %.2f ...\n', r);
    
    for s_idx = 1:length(SNR_dB)
        snr_val = SNR_lin(s_idx);
        
        % -----------------------------------------------------------------
        % 1. Theoretical PEP Calculation (Semi-Analytic)
        % -----------------------------------------------------------------
        % Effective SNR for Alamouti 2x2:
        % gamma_eff = (1/2) * SNR^(1-r) * ||H||_F^2
        snr_eff_factor = 0.5 * snr_val^(1-r);
        
        % Gamma(4,1) PDF manually implemented: f(x) = (x^3 * exp(-x)) / 6
        % This replaces 'gammapdf(x, 4, 1)' which requires Stats Toolbox
        my_gamma_pdf = @(x) (x.^3 .* exp(-x)) ./ 6;
        
        % Integration of Q-function over Gamma distribution (PEP)
        % P_e = Integral [ a*Q(sqrt(2*g*x)) * p(x) dx ]
        fun = @(x) 0.5 * erfc(sqrt(snr_eff_factor * x)/sqrt(2)) .* my_gamma_pdf(x);
        
        pep_theory = integral(fun, 0, Inf);
        
        % Store theory result (PEP approximation of SER)
        results_theory(r_idx, s_idx) = pep_theory;
        
        % -----------------------------------------------------------------
        % 2. Monte Carlo Simulation (Exact SER)
        % -----------------------------------------------------------------
        err_count = 0;
        
        % Pre-calculate scaling factors
        % Total Scaling: sqrt(0.5 * SNR^(1-r))
        scale_factor = sqrt(0.5 * snr_val^(1-r));
        
        for t = 1:num_trials
            % Channel: 2x2 Complex Gaussian
            H = (randn(2,2) + 1i*randn(2,2)) / sqrt(2);
            
            % Generate Symbol (QPSK)
            tx_idx = randi(4);
            s = qpsk_syms(tx_idx);
            
            % Alamouti Equivalent Scalar Channel
            h_fro_sq = norm(H, 'fro')^2;
            
            % Effective Noise (Unit Variance)
            n_eff = (randn + 1i*randn) / sqrt(2);
            
            % Received Signal
            y_rx = sqrt(h_fro_sq) * scale_factor * s + n_eff;
            
            % ML Detection
            [~, det_idx] = min(abs(y_rx - sqrt(h_fro_sq) * scale_factor * qpsk_syms));
            
            if det_idx ~= tx_idx
                err_count = err_count + 1;
            end
        end
        
        results_sim(r_idx, s_idx) = err_count / num_trials;
    end
end

%% =========================================================================
%% PLOTTING
%% =========================================================================
figure('Position', [100, 100, 900, 600]);

colors = ['k', 'b', 'r', 'g', 'm'];
markers = ['o', 's', 'd', '^', 'v'];

for r_idx = 1:length(r_values)
    r = r_values(r_idx);
    c = colors(r_idx);
    m = markers(r_idx);
    
    % Theoretical Slope Calculation for Legend
    % d = 4(1-r)
    slope = 4 * (1 - r);
    leg_str = sprintf('r = %.2f (Slope \\approx -%.1f)', r, slope);
    
    % Plot Simulation (Points)
    semilogy(SNR_dB, results_sim(r_idx, :), [c '-'], 'LineWidth', 1.5, ...
        'MarkerFaceColor', c, 'DisplayName', [leg_str ' Sim']);
    hold on;
    
    % Plot Theory (Dashed Lines)
    semilogy(SNR_dB, results_theory(r_idx, :), [c '--'], 'LineWidth', 2, ...
        'DisplayName', 'Theory (PEP)');
end

grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Rate (SER) / PEP');
title({'Alamouti D-M Tradeoff Visualization', 'SER Slope varies with Multiplexing Gain r'});
legend('show', 'Location', 'SouthWest');
ylim([1e-6, 1]);

% Annotation explaining the slope
% annotation('textbox', [0.6, 0.7, 0.2, 0.1], 'String', ...
%     {'Steeper Slope = Higher Diversity', ...
%      'Slope d(r) = 4(1-r)', ...
%      'r=0 -> Slope 4', ...
%      'r=1 -> Slope 0'}, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white');