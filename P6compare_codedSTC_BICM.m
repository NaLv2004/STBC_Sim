% =========================================================================
% MIMO Comparison: Alamouti (Polar Coded) vs. 5G NR Polar (SM)
% Scenario: 2x2 MIMO
%
% =========================================================================

clear; clc; close all;

%% 1. Simulation Configuration
% =========================================================================
sim_cfg.snr_db       = -4:2.0:4;          
sim_cfg.target_err   = 5000;             
sim_cfg.max_blocks   = 50000;            
sim_cfg.fading_type  = 'slow';          % 'fast' or 'slow'
sim_cfg.use_interleaver = true;        

% Polar Code Info Lengths (K)
polar_info_lens = [128]; 

% Target Spectral Efficiency (bits/channel use)
% Setting to 0.8 so Alamouti is coded (Rate ~ 0.8) and SM is coded (Rate ~ 0.4)
target_se = 0.5; 

fprintf('Running Simulation: Alamouti (Polar) vs. 5G NR Polar (SM)\n');
fprintf('Target Spectral Efficiency: %.1f bits/s/Hz\n', target_se);

%% 2. Simulations
% =========================================================================
% Storage for results
ber_alamouti_results = zeros(length(polar_info_lens), length(sim_cfg.snr_db));
ber_sm_results       = zeros(length(polar_info_lens), length(sim_cfg.snr_db));

for l_idx = 1:length(polar_info_lens)
    K = polar_info_lens(l_idx);
    
    % --- Calculate Block Lengths for Consistency ---
    % Alamouti sends 1 symbol per channel use (averaged).
    % E_alamouti = Number of encoded bits = Number of symbols (BPSK).
    E_alamouti = ceil(K / target_se);
    
    % Spatial Multiplexing sends 2 symbols per channel use.
    % To keep channel uses consistent: ChannelUses = E_alamouti.
    % E_sm = 2 * ChannelUses.
    E_sm = 2 * E_alamouti;
    
    fprintf('\n=== Info Length K = %d ===\n', K);
    fprintf('  Alamouti: E = %d (Rate %.2f)\n', E_alamouti, K/E_alamouti);
    fprintf('  SM      : E = %d (Rate %.2f)\n', E_sm, K/E_sm);
    
    % --- Run Alamouti (Polar Coded) ---
    fprintf('  Running Alamouti...\n');
    for s_idx = 1:length(sim_cfg.snr_db)
        snr_val = sim_cfg.snr_db(s_idx);
        % Speed optimization
        if s_idx > 1 && ber_alamouti_results(l_idx, s_idx-1) == 0 && ber_alamouti_results(l_idx, s_idx-1) < 1e-6
             ber_alamouti_results(l_idx, s_idx) = 0;
             continue;
        end
        ber = run_alamouti_polar(snr_val, K, E_alamouti, sim_cfg);
        ber_alamouti_results(l_idx, s_idx) = ber;
    end
    
    % --- Run Spatial Multiplexing (Polar Coded) ---
    fprintf('  Running SM...\n');
    for s_idx = 1:length(sim_cfg.snr_db)
        snr_val = sim_cfg.snr_db(s_idx);
        % Speed optimization
        if s_idx > 1 && ber_sm_results(l_idx, s_idx-1) == 0 && ber_sm_results(l_idx, s_idx-1) < 1e-6
             ber_sm_results(l_idx, s_idx) = 0;
             continue;
        end
        ber = run_5g_polar(snr_val, K, E_sm, sim_cfg);
        ber_sm_results(l_idx, s_idx) = ber;
    end
end

%% 4. Plotting
% =========================================================================
figure('Position', [100, 100, 800, 600]);
markers = {'s-', '^-', 'd-', 'v-'};
colors  = {'b', 'r', 'g', 'm'};

for l_idx = 1:length(polar_info_lens)
    K_val = polar_info_lens(l_idx);
    
    % Plot Alamouti
    semilogy(sim_cfg.snr_db, ber_alamouti_results(l_idx, :), ...
        ['--' markers{mod(l_idx-1,4)+1}(1)], 'Color', colors{mod(l_idx-1,4)+1}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Alamouti Polar (K=%d)', K_val));
    hold on;
    
    % Plot SM
    semilogy(sim_cfg.snr_db, ber_sm_results(l_idx, :), ...
        ['-' markers{mod(l_idx-1,4)+1}(1)], 'Color', colors{mod(l_idx-1,4)+1}, ...
        'LineWidth', 2, 'MarkerFaceColor', colors{mod(l_idx-1,4)+1}, ...
        'DisplayName', sprintf('SM Polar (K=%d)', K_val));
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(sprintf('2x2 MIMO: Alamouti vs SM (Polar Coded, SE=%.1f bpcu)', target_se));
legend('Location', 'SouthWest');
ylim([1e-5, 1]);

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
        H_slow = (randn(2,2)+1i*randn(2,2))/sqrt(2);
        % 1. Encode
        tx_codeword = nrPolarEncode(data, E);
        
        % 2. Interleave
        if cfg.use_interleaver
            int_idx = randperm(E).';
            tx_seq = tx_codeword(int_idx);
        else
            tx_seq = tx_codeword;
        end
        
        % 3. Modulate (BPSK)
        s_mod = 1 - 2*tx_seq;
        
        % Pad to even length for Alamouti pairing
        if mod(length(s_mod), 2) ~= 0
            s_mod = [s_mod; 0];
            was_padded = true;
        else
            was_padded = false;
        end
        
        L_syms = length(s_mod);
        num_slots = L_syms / 2;
        
        % --- CHANNEL & RECEIVER ---
        % Alamouti needs to process 2 symbols over 2 time slots
        llrs_serial = zeros(L_syms, 1);
        
        for t = 1:num_slots
            % Pair of symbols
            idx1 = 2*t-1; idx2 = 2*t;
            s1 = s_mod(idx1);
            s2 = s_mod(idx2);
            
            % Alamouti Transmission Matrix
            % Time 1: [s1; s2], Time 2: [-s2*; s1*] (conjugate irrelevant for BPSK real, but kept for rigor)
            X = [s1, -conj(s2); s2, conj(s1)];
            
            % Channel
            if strcmp(cfg.fading_type, 'slow')
                % Quasi-static over the block? Usually Alamouti assumes static over 2 slots
                % Here we generate one H per pair to be safe with 'fast' fading logic
                H = H_slow;
            else
                H = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            end
            
            N = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            Y = H * X * scale + N;
            
            % Combining
            h11=H(1,1); h12=H(1,2); h21=H(2,1); h22=H(2,2);
            y1=Y(1,1); y2=Y(2,1);   % Rx at Time 1
            y3=Y(1,2); y4=Y(2,2);   % Rx at Time 2
            
            % Standard Alamouti Combiner
            z1 = conj(h11)*y1 + conj(h21)*y2 + h12*conj(y3) + h22*conj(y4);
            z2 = conj(h12)*y1 + conj(h22)*y2 - h11*conj(y3) - h21*conj(y4);
            
            % Effective Channel Gain (Diversity Order 4 terms)
            h_sq = abs(h11)^2 + abs(h12)^2 + abs(h21)^2 + abs(h22)^2;
            
            % Soft LLR Calculation
            % z_i approx h_sq * s_i * scale + noise
            % LLR = |dist(-1)|^2 - |dist(+1)|^2  (positive -> bit 0 -> +1)
            
            % Symbol 1
            d_plus1 = abs(z1 - h_sq * (+1) * scale)^2;
            d_minus1 = abs(z1 - h_sq * (-1) * scale)^2;
            llrs_serial(idx1) = d_minus1 - d_plus1;
            
            % Symbol 2
            d_plus2 = abs(z2 - h_sq * (+1) * scale)^2;
            d_minus2 = abs(z2 - h_sq * (-1) * scale)^2;
            llrs_serial(idx2) = d_minus2 - d_plus2;
        end
        
        if was_padded, llrs_serial = llrs_serial(1:end-1); end
        
        % 5. De-Interleave
        if cfg.use_interleaver
            llrs_demod = zeros(E, 1);
            llrs_demod(int_idx) = llrs_serial;
        else
            llrs_demod = llrs_serial;
        end
        
        % 6. Decode
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
    
    % Constellation BPSK
    c_syms = 1 - 2*[0 0; 0 1; 1 0; 1 1]; 
    
    % Polar List Decoding Length
    L_list = 8;
    
    while errs < cfg.target_err && blocks < cfg.max_blocks
        % --- TRANSMITTER ---
        data = randi([0 1], K, 1);
        
        % 1. Encode
        tx_codeword = nrPolarEncode(data, E);
        
        % 2. Interleave
        if cfg.use_interleaver
            int_idx = randperm(E).';
            tx_seq = tx_codeword(int_idx);
        else
            tx_seq = tx_codeword;
        end
        
        % 3. Modulate & MIMO Map
        s_mod = 1 - 2*tx_seq; 
        
        % Padding for 2 antennas
        if mod(length(s_mod), 2) ~= 0
            s_mod = [s_mod; 0];
            was_padded = true;
        else
            was_padded = false;
        end
        
        L_syms = length(s_mod);
        X = reshape(s_mod, 2, L_syms/2);
        
        % --- CHANNEL ---
        Y = zeros(2, L_syms/2);
        scale = sqrt(snr_lin/2); 
        
        if strcmp(cfg.fading_type, 'slow')
            H = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            H_store = repmat(H, [1, 1, L_syms/2]);
        else
            H_store = (randn(2,2,L_syms/2) + 1i*randn(2,2,L_syms/2))/sqrt(2);
        end
        N_noise = (randn(2,L_syms/2) + 1i*randn(2,L_syms/2))/sqrt(2);
        
        for t = 1:L_syms/2
            Y(:,t) = H_store(:,:,t) * X(:,t) * scale + N_noise(:,t);
        end
        
        % --- RECEIVER ---
        % 4. Soft Demap (Max-Log-MAP)
        llrs_serial = zeros(L_syms, 1);
        for t = 1:L_syms/2
            y_t = Y(:,t); H_t = H_store(:,:,t);
            dists = zeros(4, 1);
            for k=1:4
                dists(k) = norm(y_t - H_t * c_syms(k,:).' * scale)^2;
            end
            % Ant 1 (Bit 1): 0->{1,2}, 1->{3,4}
            llrs_serial(2*t-1) = min(dists([3,4])) - min(dists([1,2]));
            % Ant 2 (Bit 2): 0->{1,3}, 1->{2,4}
            llrs_serial(2*t)   = min(dists([2,4])) - min(dists([1,3]));
        end
        
        if was_padded, llrs_serial = llrs_serial(1:end-1); end
        
        % 5. De-Interleave
        if cfg.use_interleaver
            llrs_demod = zeros(E, 1);
            llrs_demod(int_idx) = llrs_serial;
        else
            llrs_demod = llrs_serial;
        end
        
        % 6. Decode
        dec_data = nrPolarDecode(llrs_demod, K, E, L_list);
        
        errs = errs + sum(abs(data - double(dec_data)));
        total_bits = total_bits + K;
        blocks = blocks + 1;
    end
    ber = errs / total_bits;
end