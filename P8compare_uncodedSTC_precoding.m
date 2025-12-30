% =========================================================================
% MIMO Comparison (STC Vs Precoding with CSI)
% Scenario: 2x2 MIMO
% Scheme 1: Alamouti (No CSIT) + QPSK
% Scheme 2: MRT (Perfect CSIT) + BPSK
% Scheme 3: MRT (Perfect CSIT) + QPSK  <-- Added
% Termination: Based on minimum error count (>200) or max bits
% =========================================================================

clear; clc; close all;

fprintf('Running MIMO Comparison: Alamouti vs MRT (BPSK/QPSK)...\n');
fprintf('Simulation stops when %d bit errors are collected or %d bits sent.\n', 200, 1e7);

%% Parameters
SNR_dB = -2:2:12;            % SNR range in dB
SNR_lin = 10.^(SNR_dB./10); 

% Termination Criteria
target_errors = 1000;        % Stop after collecting this many errors
max_bits = 1e7;             % Hard limit
print_interval = 1000;      % Print status every 1000 bits

% Constellation Mappings
% QPSK (Gray Mapping): 00->1, 01->2, 10->3, 11->4
qpsk_map = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2); 
% BPSK: 0->1, 1->-1
bpsk_map = [1, -1]; 

% Result Storage
ber_alamouti_qpsk = zeros(size(SNR_dB));
ber_mrt_bpsk      = zeros(size(SNR_dB));
ber_mrt_qpsk      = zeros(size(SNR_dB));

%% Simulation Loop
for s_idx = 1:length(SNR_dB)
    snr = SNR_lin(s_idx);
    fprintf('\n--- Simulating SNR = %d dB ---\n', SNR_dB(s_idx));
    
    %% ====================================================================
    %% SCHEME 1: Alamouti (QPSK)
    %% ====================================================================
    err_bits = 0;
    total_bits = 0;
    last_print = 0;
    
    fprintf('  [Alamouti QPSK] Starting...\n');
    
    while (err_bits < target_errors) && (total_bits < max_bits)
        H = (randn(2,2) + 1i*randn(2,2)) / sqrt(2);
        
        % 4 bits -> 2 symbols
        tx_bits = randi([0, 1], 1, 4);
        
        % Map
        idx1 = 2*tx_bits(1) + tx_bits(2) + 1; 
        idx2 = 2*tx_bits(3) + tx_bits(4) + 1;
        s1 = qpsk_map(idx1);
        s2 = qpsk_map(idx2);
        
        % Channel Gain & Power Scaling
        h_pow = norm(H, 'fro')^2;
        scale = sqrt(snr / 2 * h_pow); % Power split /2
        
        % Rx
        n1 = (randn + 1i*randn)/sqrt(2);
        n2 = (randn + 1i*randn)/sqrt(2);
        y1 = scale * s1 + n1;
        y2 = scale * s2 + n2;
        
        % Detect
        [~, det_i1] = min(abs(y1 - scale * qpsk_map));
        [~, det_i2] = min(abs(y2 - scale * qpsk_map));
        
        % Demap
        bit_table = [0 0; 0 1; 1 0; 1 1];
        rx_bits = zeros(1, 4);
        rx_bits(1:2) = bit_table(det_i1, :);
        rx_bits(3:4) = bit_table(det_i2, :);
        
        % Count
        err_bits = err_bits + sum(abs(tx_bits - rx_bits));
        total_bits = total_bits + 4;
        
        if (total_bits - last_print) >= print_interval
            fprintf('  [Alamouti QPSK] Bits: %8d | Errors: %4d | BER: %.2e\n', ...
                total_bits, err_bits, err_bits/total_bits);
            last_print = total_bits;
        end
    end
    ber_alamouti_qpsk(s_idx) = err_bits / total_bits;
    fprintf('  [Alamouti QPSK] Done. Final BER: %.2e\n', ber_alamouti_qpsk(s_idx));
    
    %% ====================================================================
    %% SCHEME 2: MRT (BPSK)
    %% ====================================================================
    err_bits = 0;
    total_bits = 0;
    last_print = 0;
    
    fprintf('  [MRT BPSK     ] Starting...\n');
    
    while (err_bits < target_errors) && (total_bits < max_bits)
        H = (randn(2,2) + 1i*randn(2,2)) / sqrt(2);
        [~, S, ~] = svd(H);
        lambda_max = S(1,1);
        
        % 1 bit -> 1 symbol
        tx_bit = randi([0, 1]);
        s = bpsk_map(tx_bit + 1);
        
        % Rx (Full power on best mode)
        scale = sqrt(snr) * lambda_max;
        n = (randn + 1i*randn)/sqrt(2);
        y = scale * s + n;
        
        % Detect
        det_bit = (real(y) < 0);
        
        if det_bit ~= tx_bit
            err_bits = err_bits + 1;
        end
        total_bits = total_bits + 1;
        
        if (total_bits - last_print) >= print_interval
            fprintf('  [MRT BPSK     ] Bits: %8d | Errors: %4d | BER: %.2e\n', ...
                total_bits, err_bits, err_bits/total_bits);
            last_print = total_bits;
        end
    end
    ber_mrt_bpsk(s_idx) = err_bits / total_bits;
    fprintf('  [MRT BPSK     ] Done. Final BER: %.2e\n', ber_mrt_bpsk(s_idx));
    
    %% ====================================================================
    %% SCHEME 3: MRT (QPSK)  <-- NEW ADDITION
    %% ====================================================================
    err_bits = 0;
    total_bits = 0;
    last_print = 0;
    
    fprintf('  [MRT QPSK     ] Starting...\n');
    
    while (err_bits < target_errors) && (total_bits < max_bits)
        H = (randn(2,2) + 1i*randn(2,2)) / sqrt(2);
        [~, S, ~] = svd(H);
        lambda_max = S(1,1);
        
        % 2 bits -> 1 symbol
        tx_bits = randi([0, 1], 1, 2);
        
        % Map
        idx = 2*tx_bits(1) + tx_bits(2) + 1;
        s = qpsk_map(idx);
        
        % Rx (Full power on best mode)
        scale = sqrt(snr) * lambda_max;
        n = (randn + 1i*randn)/sqrt(2);
        y = scale * s + n;
        
        % Detect
        [~, det_idx] = min(abs(y - scale * qpsk_map));
        
        % Demap
        bit_table = [0 0; 0 1; 1 0; 1 1];
        rx_bits = bit_table(det_idx, :);
        
        % Count
        err_bits = err_bits + sum(abs(tx_bits - rx_bits));
        total_bits = total_bits + 2;
        
        if (total_bits - last_print) >= print_interval
            fprintf('  [MRT QPSK     ] Bits: %8d | Errors: %4d | BER: %.2e\n', ...
                total_bits, err_bits, err_bits/total_bits);
            last_print = total_bits;
        end
    end
    ber_mrt_qpsk(s_idx) = err_bits / total_bits;
    fprintf('  [MRT QPSK     ] Done. Final BER: %.2e\n', ber_mrt_qpsk(s_idx));
end

%% Plotting
figure('Position', [200, 200, 800, 600]);

% 1. Alamouti QPSK (Blue Square)
semilogy(SNR_dB, ber_alamouti_qpsk, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;

% 2. MRT BPSK (Red Circle)
semilogy(SNR_dB, ber_mrt_bpsk, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');

% 3. MRT QPSK (Green Triangle)
semilogy(SNR_dB, ber_mrt_qpsk, 'g-^', 'LineWidth', 2, 'MarkerFaceColor', 'g');

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title({'2x2 MIMO Performance Comparison', 'Effect of CSI and Modulation'});
legend('Alamouti (QPSK, No CSI)', ...
       'MRT (BPSK, Perfect CSI)', ...
       'MRT (QPSK, Perfect CSI)', ...
       'Location', 'SouthWest');
ylim([1e-6, 1]);

fprintf('\nAll Simulations Complete.\n');