% =========================================================================
% MIMO Comparison: Uncoded Alamouti vs. 5G NR LDPC (Fixed)
% Scenario: 2x2 MIMO, Rate 1 bit/s/Hz
% 
% =========================================================================

clear; clc; close all;

%% 1. Simulation Configuration
% =========================================================================
sim_cfg.snr_db       = -4:0.5:0.5;          
sim_cfg.target_err   = 5000;             
sim_cfg.max_blocks   = 20000;            
sim_cfg.fading_type  = 'fast';          % 'fast' or 'slow'
sim_cfg.use_interleaver = true;        

% 5G Standard K values (must be multiples of 10 or 22 depending on BG)
% We use a helper to snap these to valid values automatically.
ldpc_info_lens_input = [128, 256, 528, 1056]; 

fprintf('Running Simulation: Alamouti vs. 5G NR LDPC (Manual Puncturing Fix)\n');

%% 2. Run Alamouti Baseline
% =========================================================================
fprintf('\n--- Running Scheme 1: Alamouti (Uncoded) ---\n');
ber_alamouti = zeros(size(sim_cfg.snr_db));
for i = 1:length(sim_cfg.snr_db)
    ber_alamouti(i) = run_alamouti(sim_cfg.snr_db(i), sim_cfg.target_err, 1e5);
    fprintf('  SNR %2d dB | BER: %.2e\n', sim_cfg.snr_db(i), ber_alamouti(i));
end

%% 3. Run 5G LDPC Schemes
% =========================================================================
ber_ldpc_results = zeros(length(ldpc_info_lens_input), length(sim_cfg.snr_db));

for l_idx = 1:length(ldpc_info_lens_input)
    K_req = ldpc_info_lens_input(l_idx);
    R = 0.5;
    
    % Get valid K, BaseGraph, and Zc (Lifting Size)
    [K, bgn, Zc] = get_5g_params(K_req, R);
    
    if K ~= K_req
        fprintf('\n[Info] Adjusted K from %d to %d (BG%d, Zc=%d)\n', K_req, K, bgn, Zc);
    end
    
    fprintf('--- Running Scheme 2: SM + 5G LDPC (K=%d, Rate=0.5, BG%d) ---\n', K, bgn);
    
    for s_idx = 1:length(sim_cfg.snr_db)
        snr_val = sim_cfg.snr_db(s_idx);
        
        % Speed optimization
        if s_idx > 1 && ber_ldpc_results(l_idx, s_idx-1) == 0
            ber_ldpc_results(l_idx, s_idx) = 0;
            fprintf('  SNR %2d dB | BER: <Skipped>\n', snr_val);
            continue;
        end
        
        ber = run_5g_ldpc_manual(snr_val, K, R, bgn, Zc, sim_cfg);
        ber_ldpc_results(l_idx, s_idx) = ber;
        fprintf('  SNR %2d dB | BER: %.2e\n', snr_val, ber);
    end
end

%% 4. Plotting
% =========================================================================
figure('Position', [100, 100, 800, 600]);
markers = {'s-', '^-', 'd-', 'v-'};
colors  = {'b', 'r', 'g', 'm'};

semilogy(sim_cfg.snr_db, ber_alamouti, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Alamouti (Uncoded)');
hold on;

for l_idx = 1:length(ldpc_info_lens_input)
    [K_used, ~] = get_5g_params(ldpc_info_lens_input(l_idx), 0.5);
    legend_str = sprintf('SM + 5G LDPC (K=%d, R=0.5)', K_used);
    semilogy(sim_cfg.snr_db, ber_ldpc_results(l_idx, :), ...
        markers{mod(l_idx-1,4)+1}, 'Color', colors{mod(l_idx-1,4)+1}, ...
        'LineWidth', 2, 'MarkerFaceColor', colors{mod(l_idx-1,4)+1}, ...
        'DisplayName', legend_str);
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('2x2 MIMO: Alamouti vs 5G LDPC (Corrected)');
legend('Location', 'SouthWest');
ylim([1e-6, 1]);

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function [K_valid, bgn, Zc] = get_5g_params(K_req, R)
    % 1. Determine Base Graph (BG) per 3GPP TS 38.212
    if K_req <= 292 || (K_req <= 3824 && R <= 0.67) || R <= 0.25
        bgn = 2;
        kb = 10; % K = 10 * Zc
    else
        bgn = 1;
        kb = 22; % K = 22 * Zc
    end
    
    % 2. Find valid Lifting Size (Zc)
    % 5G NR defined Zc values
    Zc_all = [2:16, 18:2:32, 36:4:64, 72:8:128, 144:16:256, 288:32:384];
    valid_Ks = kb * Zc_all;
    
    % Snap to nearest valid K
    [~, idx] = min(abs(valid_Ks - K_req));
    Zc = Zc_all(idx);
    K_valid = valid_Ks(idx);
end

function ber = run_alamouti(snr_db, target_err, max_bits)
    snr_lin = 10^(snr_db/10);
    errs = 0; bits = 0;
    scale = sqrt(snr_lin/2); 
    
    while errs < target_err && bits < max_bits
        H = (randn(2,2)+1i*randn(2,2))/sqrt(2);
        tx_bits = randi([0 1], 1, 2);
        s = 1 - 2*tx_bits; 
        X = [s(1), -conj(s(2)); s(2), conj(s(1))];
        N = (randn(2,2)+1i*randn(2,2))/sqrt(2);
        Y = H * X * scale + N;
        h11=H(1,1); h12=H(1,2); h21=H(2,1); h22=H(2,2);
        y1=Y(1,1); y2=Y(2,1); y3=Y(1,2); y4=Y(2,2);
        z1 = conj(h11)*y1 + conj(h21)*y2 + h12*conj(y3) + h22*conj(y4);
        z2 = conj(h12)*y1 + conj(h22)*y2 - h11*conj(y3) - h21*conj(y4);
        det_bits = [real(z1)<0, real(z2)<0];
        errs = errs + sum(abs(tx_bits - det_bits));
        bits = bits + 2;
    end
    ber = errs/bits;
end

function ber = run_5g_ldpc_manual(snr_db, K, R, bgn, Zc, cfg)
    snr_lin = 10^(snr_db/10);
    errs = 0; total_bits = 0; blocks = 0;
    
    % Target transmission length
    E = ceil(K / R); 
    
    % Constellation BPSK
    c_syms = 1 - 2*[0 0; 0 1; 1 0; 1 1]; 
    
    while errs < cfg.target_err && blocks < cfg.max_blocks
        % --- TRANSMITTER ---
        data = randi([0 1], K, 1);
        
        % 1. Encode: Get Full Circular Buffer (Length N)
        enc_bits_full = nrLDPCEncode(data, bgn);
        N_cw = length(enc_bits_full);
        
        % 2. Manual Rate Matching (Crucial Step)
        % 5G Spec: The first 2*Zc bits are ALWAYS punctured (Systematic).
        % We must start transmission from index 2*Zc + 1.
        start_idx = 2 * Zc + 1;
        
        % Check if we wrap around (for R=0.5 we usually don't, but safe to check)
        indices = mod((start_idx-1 : start_idx+E-2), N_cw) + 1;
        tx_codeword = enc_bits_full(indices);
        
        % 3. Interleave
        if cfg.use_interleaver
            int_idx = randperm(E).';
            tx_seq = tx_codeword(int_idx);
        else
            tx_seq = tx_codeword;
        end
        
        % 4. Modulate & MIMO Map
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
        % 5. Soft Demap (Max-Log-MAP)
        % LLR = ln(P(0)/P(1)). Positive -> 0.
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
        
        % 6. De-Interleave
        if cfg.use_interleaver
            llrs_demod = zeros(E, 1);
            llrs_demod(int_idx) = llrs_serial;
        else
            llrs_demod = llrs_serial;
        end
        
        % 7. Manual Rate Recover (The Fix)
        % Create a buffer of zeros (LLR=0 means uncertainty/punctured)
        llrs_decoder_input = zeros(N_cw, 1);
        
        % Place the received LLRs into the exact positions we took them from
        llrs_decoder_input(indices) = llrs_demod;
        
        % Note: The first 2*Zc bits remain 0 (Correctly handled as punctured)
        
        % 8. Decode
        % DecisionType='hard' returns bits.
        dec_bits = nrLDPCDecode(llrs_decoder_input, bgn, 25, ...
            'DecisionType', 'hard', 'Algorithm', 'Layered Belief Propagation'); 
        
        dec_data = dec_bits(1:K);
        
        errs = errs + sum(abs(data - double(dec_data)));
        total_bits = total_bits + K;
        blocks = blocks + 1;
    end
    ber = errs / total_bits;
end