% =========================================================================
% MIMO Comparison: Alamouti (Polar Coded) vs. 5G NR Polar (SM)
% Scenario: 2x2 MIMO with Variable Fading Speed
% Note: Enforced Block Fading alignment for Alamouti pairs.
% =========================================================================

clear; clc; close all;

%% 1. Simulation Configuration
% =========================================================================
sim_cfg.snr_db       = -8:1:4;          
sim_cfg.target_err   = 1000;             
sim_cfg.max_blocks   = 10000;            
sim_cfg.use_interleaver = true;        

% Polar Code Info Lengths (K)
polar_info_lens = [64,128]; 

% Fading Parameters (K denotes channel invariant for K channel uses/time slots)
% 2   = Channel constant for exactly 1 Alamouti pair (Fastest supported fading)
% 10  = Channel constant for 5 Alamouti pairs
% 100 = Slow fading
fading_param_list = [2, 16, 64]; 

% Target Spectral Efficiency (bits/channel use)
target_se = 0.5; 

fprintf('Running Simulation: Alamouti (Polar) vs. 5G NR Polar (SM)\n');
fprintf('Target Spectral Efficiency: %.1f bits/cu\n', target_se);
fprintf('Fading Parameters (Coherence Time in Slots): %s\n', num2str(fading_param_list));

%% 2. Simulations
% =========================================================================
results_idx = 0;
results = struct([]);

for f_idx = 1:length(fading_param_list)
    current_T_coh = fading_param_list(f_idx);
    sim_cfg.T_coh = current_T_coh; 
    
    fprintf('\n==================================================\n');
    fprintf('Running Fading Scenario: Channel Coherence = %d Time Slots\n', current_T_coh);
    fprintf('==================================================\n');

    for l_idx = 1:length(polar_info_lens)
        K = polar_info_lens(l_idx);
        
        % --- Block Length Calculation ---
        % Alamouti: Rate = K / E_alamouti (bits per symbol)
        % Since Alamouti takes 2 slots for 2 symbols, Rate per Channel Use = Rate.
        E_alamouti = ceil(K / target_se);
        
        % SM: Rate = K / E_sm (bits per symbol).
        % SM takes 1 slot for 2 symbols. To match SE, we need consistent total slots.
        % E_sm = 2 * E_alamouti.
        E_sm = 2 * E_alamouti;
        
        fprintf('\n[Info Length K = %d]\n', K);
        
        % --- Run Alamouti ---
        fprintf('  > Scheme: Alamouti Polar (Rate %.2f)\n', K/E_alamouti);
        ber_temp = zeros(size(sim_cfg.snr_db));
        for s_idx = 1:length(sim_cfg.snr_db)
            snr_val = sim_cfg.snr_db(s_idx);
            
            % Speed optimization
            if s_idx > 1 && ber_temp(s_idx-1) == 0
                 ber_temp(s_idx) = 0;
                 fprintf('    SNR %5.1f dB | BER: < 1e-6 (Skipped)\n', snr_val);
                 continue;
            end
            
            ber = run_alamouti_polar(snr_val, K, E_alamouti, sim_cfg);
            ber_temp(s_idx) = ber;
            fprintf('    SNR %5.1f dB | BER: %.2e\n', snr_val, ber);
        end
        
        % Store Alamouti
        results_idx = results_idx + 1;
        results(results_idx).type = 'Alamouti';
        results(results_idx).K = K;
        results(results_idx).T_coh = current_T_coh;
        results(results_idx).snr = sim_cfg.snr_db;
        results(results_idx).ber = ber_temp;
        results(results_idx).color = 'k'; 
        results(results_idx).marker = 's';
        
        % --- Run SM ---
        fprintf('  > Scheme: SM Polar (Rate %.2f)\n', K/E_sm);
        ber_temp = zeros(size(sim_cfg.snr_db));
        for s_idx = 1:length(sim_cfg.snr_db)
            snr_val = sim_cfg.snr_db(s_idx);
            
            if s_idx > 1 && ber_temp(s_idx-1) == 0
                 ber_temp(s_idx) = 0;
                 fprintf('    SNR %5.1f dB | BER: < 1e-6 (Skipped)\n', snr_val);
                 continue;
            end
            
            ber = run_5g_polar(snr_val, K, E_sm, sim_cfg);
            ber_temp(s_idx) = ber;
            fprintf('    SNR %5.1f dB | BER: %.2e\n', snr_val, ber);
        end
        
        % Store SM
        results_idx = results_idx + 1;
        results(results_idx).type = 'SM';
        results(results_idx).K = K;
        results(results_idx).T_coh = current_T_coh;
        results(results_idx).snr = sim_cfg.snr_db;
        results(results_idx).ber = ber_temp;
        results(results_idx).color = 'b';
        results(results_idx).marker = 'o';
    end
end

%% 4. Plotting
% =========================================================================
figure('Position', [100, 100, 900, 700]);

% Line styles for different fading speeds
line_styles = {'-', '--', '-.', ':'}; 
colors_ala  = [0 0 0; 0.4 0.4 0.4; 0.7 0.7 0.7]; % Black gradients
colors_sm   = [0 0 1; 0.3 0.6 1; 0.6 0.8 1];     % Blue gradients

hold on;
legend_entries = {};
unique_T = sort(unique([results.T_coh]));

for i = 1:length(results)
    res = results(i);
    
    % Determine style based on T_coh index
    t_idx = find(unique_T == res.T_coh);
    l_style = line_styles{mod(t_idx-1, length(line_styles)) + 1};
    
    % Determine color
    if strcmp(res.type, 'Alamouti')
        c_val = colors_ala(mod(t_idx-1, size(colors_ala,1))+1, :);
    else
        c_val = colors_sm(mod(t_idx-1, size(colors_sm,1))+1, :);
    end
    
    semilogy(res.snr, res.ber, ...
        'LineStyle', l_style, ...
        'Color', c_val, ...
        'Marker', res.marker, ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', c_val);
    
    legend_entries{end+1} = sprintf('%s (T_{coh}=%d)', res.type, res.T_coh);
end

grid on; grid minor;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(sprintf('2x2 MIMO Comparison: Target SE = %.1f bpcu', target_se), 'FontSize', 14);
legend(legend_entries, 'Location', 'SouthWest', 'FontSize', 10);
ylim([1e-6, 1]);

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function ber = run_alamouti_polar(snr_db, K, E, cfg)
    snr_lin = 10^(snr_db/10);
    errs = 0; total_bits = 0; blocks = 0;
    L_list = 8;
    scale = sqrt(snr_lin/2); 
    
    while errs < cfg.target_err && blocks < cfg.max_blocks
        % --- TRANSMITTER ---
        data = randi([0 1], K, 1);
        tx_codeword = nrPolarEncode(data, E);
        
        if cfg.use_interleaver
            int_idx = randperm(E).';
            tx_seq = tx_codeword(int_idx);
        else
            tx_seq = tx_codeword;
        end
        
        s_mod = 1 - 2*tx_seq; % BPSK
        
        if mod(length(s_mod), 2) ~= 0
            s_mod = [s_mod; 0];
            was_padded = true;
        else
            was_padded = false;
        end
        
        L_syms = length(s_mod);
        num_slots = L_syms / 2;
        
        % --- CHANNEL GENERATION ---
        T_coh = cfg.T_coh;
        
        % Generate enough fading blocks for the whole frame
        num_fading_blocks = ceil(L_syms / T_coh);
        H_blocks = (randn(2,2,num_fading_blocks) + 1i*randn(2,2,num_fading_blocks))/sqrt(2);
        
        llrs_serial = zeros(L_syms, 1);
        
        for t = 1:num_slots
            % Pair indices
            idx1 = 2*t-1; % Time Slot T
            idx2 = 2*t;   % Time Slot T+1
            
            s1 = s_mod(idx1);
            s2 = s_mod(idx2);
            
            % --- CRITICAL CHANGE FOR ALAMOUTI STABILITY ---
            % We calculate the block index based on the FIRST slot of the pair.
            % We enforce that both slots in this pair use the SAME channel matrix.
            % This guarantees H(t) == H(t+1) for the pair, preserving orthogonality
            % as long as T_coh is even (or at least >= 2).
            block_idx = ceil(idx1 / T_coh);
            
            H_pair = H_blocks(:,:,block_idx); 
            
            % Transmission (Physical Channel Application)
            % Time 1: H_pair * [s1; s2]
            x1 = [s1; s2];
            y1_vec = H_pair * x1 * scale + (randn(2,1)+1i*randn(2,1))/sqrt(2);
            
            % Time 2: H_pair * [-s2*; s1*]
            x2 = [-conj(s2); conj(s1)];
            y2_vec = H_pair * x2 * scale + (randn(2,1)+1i*randn(2,1))/sqrt(2);
            
            % --- RECEIVER COMBINING ---
            % Since we forced H1 == H2 == H_pair, standard Alamouti works perfectly here.
            
            h11=H_pair(1,1); h12=H_pair(1,2); h21=H_pair(2,1); h22=H_pair(2,2);
            r1=y1_vec(1); r2=y1_vec(2);   
            r3=y2_vec(1); r4=y2_vec(2);   
            
            % Standard Linear Combiner
            z1 = conj(h11)*r1 + conj(h21)*r2 + h12*conj(r3) + h22*conj(r4);
            z2 = conj(h12)*r1 + conj(h22)*r2 - h11*conj(r3) - h21*conj(r4);
            
            h_sq = abs(h11)^2 + abs(h12)^2 + abs(h21)^2 + abs(h22)^2;
            
            % LLRs
            d_plus1 = abs(z1 - h_sq * (+1) * scale)^2;
            d_minus1 = abs(z1 - h_sq * (-1) * scale)^2;
            llrs_serial(idx1) = d_minus1 - d_plus1;
            
            d_plus2 = abs(z2 - h_sq * (+1) * scale)^2;
            d_minus2 = abs(z2 - h_sq * (-1) * scale)^2;
            llrs_serial(idx2) = d_minus2 - d_plus2;
        end
        
        if was_padded, llrs_serial = llrs_serial(1:end-1); end
        
        if cfg.use_interleaver
            llrs_demod = zeros(E, 1);
            llrs_demod(int_idx) = llrs_serial;
        else
            llrs_demod = llrs_serial;
        end
        
        dec_data = nrPolarDecode(llrs_demod, K, E, L_list);
        errs = errs + sum(abs(data - double(dec_data)));
        total_bits = total_bits + K;
        blocks = blocks + 1;
    end
    ber = errs/total_bits;
end

function ber = run_5g_polar(snr_db, K, E, cfg)
    snr_lin = 10^(snr_db/10);
    errs = 0; total_bits = 0; blocks = 0;
    c_syms = 1 - 2*[0 0; 0 1; 1 0; 1 1]; 
    L_list = 8;
    
    while errs < cfg.target_err && blocks < cfg.max_blocks
        data = randi([0 1], K, 1);
        tx_codeword = nrPolarEncode(data, E);
        
        if cfg.use_interleaver
            int_idx = randperm(E).';
            tx_seq = tx_codeword(int_idx);
        else
            tx_seq = tx_codeword;
        end
        
        s_mod = 1 - 2*tx_seq; 
        if mod(length(s_mod), 2) ~= 0, s_mod = [s_mod; 0]; was_padded = true; else, was_padded = false; end
        
        L_syms = length(s_mod);
        X = reshape(s_mod, 2, L_syms/2);
        
        % --- CHANNEL GENERATION ---
        T_coh = cfg.T_coh;
        num_uses = L_syms/2; 
        num_fading_blocks = ceil(num_uses / T_coh);
        H_blocks = (randn(2,2,num_fading_blocks) + 1i*randn(2,2,num_fading_blocks))/sqrt(2);
        
        Y = zeros(2, num_uses);
        scale = sqrt(snr_lin/2); 
        
        % Apply Channel (Time Slot based)
        for t = 1:num_uses
            block_idx = ceil(t / T_coh);
            H_curr = H_blocks(:,:,block_idx);
            Y(:,t) = H_curr * X(:,t) * scale + (randn(2,1) + 1i*randn(2,1))/sqrt(2);
        end
        
        % --- RECEIVER ---
        llrs_serial = zeros(L_syms, 1);
        for t = 1:num_uses
            y_t = Y(:,t); 
            block_idx = ceil(t / T_coh);
            H_t = H_blocks(:,:,block_idx); 
            
            dists = zeros(4, 1);
            for k=1:4
                dists(k) = norm(y_t - H_t * c_syms(k,:).' * scale)^2;
            end
            llrs_serial(2*t-1) = min(dists([3,4])) - min(dists([1,2]));
            llrs_serial(2*t)   = min(dists([2,4])) - min(dists([1,3]));
        end
        
        if was_padded, llrs_serial = llrs_serial(1:end-1); end
        
        if cfg.use_interleaver
            llrs_demod = zeros(E, 1);
            llrs_demod(int_idx) = llrs_serial;
        else
            llrs_demod = llrs_serial;
        end
        
        dec_data = nrPolarDecode(llrs_demod, K, E, L_list);
        errs = errs + sum(abs(data - double(dec_data)));
        total_bits = total_bits + K;
        blocks = blocks + 1;
    end
    ber = errs / total_bits;
end