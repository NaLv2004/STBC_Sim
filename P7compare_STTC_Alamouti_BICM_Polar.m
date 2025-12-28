% =========================================================================
% MIMO Comparison: Alamouti vs. 5G NR Polar vs. STTC
% Scenario: 2x2 MIMO, Rate 1 bit/s/Hz, BPSK
%
% =========================================================================

clear; clc; close all;

%% 1. Simulation Configuration
% =========================================================================
sim_cfg.snr_db       = -4:1.0:6;          % Extended range to see STTC slope
sim_cfg.target_err   = 5000;             
sim_cfg.max_blocks   = 50000;            
sim_cfg.fading_type  = 'fast';          % 'fast' implies Time Diversity available
sim_cfg.use_interleaver = true;        

% Polar Code Info Lengths (K)
polar_info_lens = [64]; 

fprintf('Running Simulation: Alamouti vs. 5G Polar vs. STTC\n');
fprintf('Channel: 2x2 MIMO, Rayleigh Fading (%s)\n', sim_cfg.fading_type);

%% 2. Run Alamouti Baseline
% =========================================================================
fprintf('\n--- Scheme 1: Alamouti (Uncoded) ---\n');
ber_alamouti = zeros(size(sim_cfg.snr_db));
for i = 1:length(sim_cfg.snr_db)
    ber_alamouti(i) = run_alamouti(sim_cfg.snr_db(i), sim_cfg.target_err, 1e5, sim_cfg.fading_type);
    fprintf('  SNR %2d dB | BER: %.2e\n', sim_cfg.snr_db(i), ber_alamouti(i));
end

%% 3. Run 5G Polar Schemes
% =========================================================================
ber_polar_results = zeros(length(polar_info_lens), length(sim_cfg.snr_db));

for l_idx = 1:length(polar_info_lens)
    K = polar_info_lens(l_idx);
    R = 0.5; % Rate 0.5 per stream * 2 streams = Total Rate 1 (SM)
    E = ceil(K / R); 
    
    fprintf('\n--- Scheme 2: SM + 5G Polar (K=%d, Rate=1 total) ---\n', K);
    
    for s_idx = 1:length(sim_cfg.snr_db)
        snr_val = sim_cfg.snr_db(s_idx);
        if s_idx > 1 && ber_polar_results(l_idx, s_idx-1) < 1e-6
            ber_polar_results(l_idx, s_idx) = 1e-7; continue; 
        end
        ber = run_5g_polar(snr_val, K, E, sim_cfg);
        ber_polar_results(l_idx, s_idx) = ber;
        fprintf('  SNR %2d dB | BER: %.2e\n', snr_val, ber);
    end
end

%% 4. Run STTC (Space-Time Trellis Code)
% =========================================================================
fprintf('\n--- Scheme 3: 4-State STTC (BPSK) ---\n');
ber_sttc = zeros(size(sim_cfg.snr_db));
for i = 1:length(sim_cfg.snr_db)
    % STTC is computationally heavier, reduce target errors slightly for speed if needed
    ber_sttc(i) = run_sttc(sim_cfg.snr_db(i), sim_cfg.target_err, sim_cfg.max_blocks, sim_cfg.fading_type);
    fprintf('  SNR %2d dB | BER: %.2e\n', sim_cfg.snr_db(i), ber_sttc(i));
end

%% 5. Plotting
% =========================================================================
figure('Position', [100, 100, 900, 600]);
markers = {'s-', '^-', 'd-', 'v-'};

% Plot Alamouti
semilogy(sim_cfg.snr_db, ber_alamouti, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Alamouti (STBC)');
hold on;

% Plot STTC
semilogy(sim_cfg.snr_db, ber_sttc, 'b-d', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', '4-State STTC');

% Plot Polar
for l_idx = 1:length(polar_info_lens)
    K_val = polar_info_lens(l_idx);
    legend_str = sprintf('5G Polar (FEC, K=%d)', K_val);
    semilogy(sim_cfg.snr_db, ber_polar_results(l_idx, :), ...
        'r-^', 'LineWidth', 2, 'MarkerFaceColor', 'r', ...
        'DisplayName', legend_str);
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(['2x2 MIMO: STBC vs STTC vs Polar (' sim_cfg.fading_type ' fading)']);
legend('Location', 'SouthWest');
ylim([1e-6, 1]);

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

%% --- 1. Alamouti ---
function ber = run_alamouti(snr_db, target_err, max_bits, fading_type)
    snr_lin = 10^(snr_db/10);
    errs = 0; bits = 0;
    scale = sqrt(snr_lin/2); % Normalize Tx energy
    
    while errs < target_err && bits < max_bits
        tx_bits = randi([0 1], 1, 2);
        s = 1 - 2*tx_bits; % BPSK
        
        % Alamouti Coding
        X = [s(1), -conj(s(2)); s(2), conj(s(1))];
        
        % Channel
        if strcmp(fading_type, 'fast')
            % Channel changes every symbol (col of X)
            H1 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            H2 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            % Received signal
            N = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            Y = zeros(2,2);
            Y(:,1) = H1 * X(:,1) * scale + N(:,1);
            Y(:,2) = H2 * X(:,2) * scale + N(:,2);
            
            % For fast fading Alamouti, perfect orthogonality is lost unless 
            % sophisticated combining is used. Standard Alamouti assumes H const 
            % over 2 symbols. Here we approximate H_eff = (H1+H2)/2 for decoding
            % OR strictly follow standard Alamouti which fails in very fast fading.
            % To be fair to STBC, we usually assume Quasi-Static (block of 2).
            % Let's assume H is constant for the pair even in 'fast' mode for STBC fairness.
            H = H1; 
            Y = H * X * scale + N;
        else
            H = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            N = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            Y = H * X * scale + N;
        end

        % Combining
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

%% --- 2. 5G NR Polar ---
function ber = run_5g_polar(snr_db, K, E, cfg)
    snr_lin = 10^(snr_db/10);
    errs = 0; total_bits = 0; blocks = 0;
    c_syms = 1 - 2*[0 0; 0 1; 1 0; 1 1]; % QPSK mapping for distance calc
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
        X = reshape(s_mod, 2, L_syms/2); % Spatial Multiplexing
        
        scale = sqrt(snr_lin/2); 
        Y = zeros(2, L_syms/2);
        
        % Channel Generation
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
        
        % Soft Demap (Max-Log-MAP)
        llrs_serial = zeros(L_syms, 1);
        for t = 1:L_syms/2
            y_t = Y(:,t); H_t = H_store(:,:,t);
            dists = zeros(4, 1);
            % Enumerate all 4 combinations of 2 BPSK antennas
            % 00, 01, 10, 11
            for k=1:4
                dists(k) = norm(y_t - H_t * c_syms(k,:).' * scale)^2;
            end
            llrs_serial(2*t-1) = min(dists([3,4])) - min(dists([1,2]));
            llrs_serial(2*t)   = min(dists([2,4])) - min(dists([1,3]));
        end
        
        if was_padded, llrs_serial = llrs_serial(1:end-1); end
        if cfg.use_interleaver, llrs_demod = zeros(E, 1); llrs_demod(int_idx) = llrs_serial; else, llrs_demod = llrs_serial; end
        
        dec_data = nrPolarDecode(llrs_demod, K, E, L_list);
        errs = errs + sum(abs(data - double(dec_data)));
        total_bits = total_bits + K;
        blocks = blocks + 1;
    end
    ber = errs / total_bits;
end

%% --- 3. STTC (Space-Time Trellis Code) ---
function ber = run_sttc(snr_db, target_err, max_blocks, fading_type)
    % --- 修改点：使用强卷积码 (K=7, 64状态) ---
    % 生成多项式 [133 171] 是 K=7 时的最佳距离谱码
    % 状态数 = 2^(K-1) = 64
    constraint_length = 7;
    polys = [133 171]; 
    trellis = poly2trellis(constraint_length, polys);
    
    num_states = trellis.numStates; % 变为 64
    
    % STTC Rate 1 bit/s/Hz (1 input bit -> 2 output bits -> 2 Antennas BPSK)
    
    snr_lin = 10^(snr_db/10);
    scale = sqrt(snr_lin/2); 
    
    L_frame = 100; 
    
    errs = 0; total_bits = 0; blocks = 0;
    
    while errs < target_err && blocks < max_blocks
        % 1. Generate Bits
        % 为了让 Viterbi 完美收敛，通常需要加尾比特 (Tail bits)
        % 这里为了保持块长一致性，我们在数据后补 K-1 个 0
        num_data_bits = L_frame - (constraint_length - 1);
        data_bits = randi([0 1], num_data_bits, 1);
        tail_bits = zeros(constraint_length - 1, 1);
        full_bits = [data_bits; tail_bits];
        
        % 2. Encode
        enc_bits = convenc(full_bits, trellis);
        tx_syms_idx = reshape(enc_bits, 2, []).'; % [L_frame x 2]
        
        % 3. Map to BPSK
        X_frame = 1 - 2 * tx_syms_idx; 
        X_frame = X_frame.';           
        
        % 4. Channel
        if strcmp(fading_type, 'fast')
            H_store = (randn(2,2,L_frame) + 1i*randn(2,2,L_frame))/sqrt(2);
        else
            H_blk = (randn(2,2)+1i*randn(2,2))/sqrt(2);
            H_store = repmat(H_blk, [1, 1, L_frame]);
        end
        N_noise = (randn(2,L_frame) + 1i*randn(2,L_frame))/sqrt(2);
        
        Y = zeros(2, L_frame);
        for t = 1:L_frame
            Y(:,t) = H_store(:,:,t) * X_frame(:,t) * scale + N_noise(:,t);
        end
        
        % 5. MIMO Viterbi Decoder (Soft Input)
        path_metrics = inf(num_states, 1);
        path_metrics(1) = 0; 
        
        survivor_paths = zeros(num_states, L_frame);
        
        % 预计算输出映射表以加速 (BPSK映射)
        % trellis.outputs 是十进制 0-3
        % 预先转为符号值: row index = output_dec + 1
        % Map: 0->[1 1], 1->[1 -1], 2->[-1 1], 3->[-1 -1]
        % 注意：convenc 的输出高位对应多项式列表的第一个
        % 我们建立一个查找表 output_map(output_dec+1, :)
        output_map = zeros(4, 2);
        for out_dec = 0:3
            b1 = floor(out_dec/2);
            b2 = mod(out_dec, 2);
            output_map(out_dec+1, :) = [1 - 2*b1, 1 - 2*b2];
        end
        
        for t = 1:L_frame
            y_t = Y(:,t);
            H_t = H_store(:,:,t);
            new_path_metrics = inf(num_states, 1);
            
            % 优化：只遍历当前活跃的状态（Metric不为inf的）
            % 在开始阶段能显著加速
            active_states = find(path_metrics < inf);
            
            for i = 1:length(active_states)
                s = active_states(i);
                
                % Try input 0 and 1
                for input_bit = 0:1
                    next_s = trellis.nextStates(s, input_bit+1) + 1; 
                    out_dec = trellis.outputs(s, input_bit+1);       
                    
                    x_expected = output_map(out_dec+1, :).';
                    
                    bm = norm(y_t - H_t * x_expected * scale)^2;
                    
                    if path_metrics(s) + bm < new_path_metrics(next_s)
                        new_path_metrics(next_s) = path_metrics(s) + bm;
                        survivor_paths(next_s, t) = s; 
                    end
                end
            end
            path_metrics = new_path_metrics;
        end
        
        % 6. Traceback
        % 因为加了尾比特，理论上最后一定回到状态 1 (全0状态)
        % 但为了保险，还是取最小度量
        [~, best_final_state] = min(path_metrics);
        curr_s = best_final_state;
        decoded_bits_full = zeros(L_frame, 1);
        
        for t = L_frame:-1:1
            prev_s = survivor_paths(curr_s, t);
            if (trellis.nextStates(prev_s, 1) + 1) == curr_s
                decoded_bits_full(t) = 0;
            else
                decoded_bits_full(t) = 1;
            end
            curr_s = prev_s;
        end
        
        % 去除尾比特进行误码统计
        dec_data = decoded_bits_full(1:num_data_bits);
        
        errs = errs + sum(abs(data_bits - dec_data));
        total_bits = total_bits + num_data_bits;
        blocks = blocks + 1;
    end
    ber = errs / total_bits;
end
% function ber = run_sttc(snr_db, target_err, max_blocks, fading_type)
%     % 4-State STTC for 2 Tx antennas using BPSK
%     % Equivalent to Rate 1/2 Convolutional Code [5, 7] mapped to Ant1, Ant2
%     % Rate = 1 bit/s/Hz
% 
%     snr_lin = 10^(snr_db/10);
%     scale = sqrt(snr_lin/2); % Total Tx energy = 2 (1 per ant). Scale adjusts SNR.
% 
%     % Trellis Definition: Generator [5 7] octal -> [101; 111] binary
%     % Constraint Length 3 -> 4 States (00, 10, 01, 11) -> 0, 1, 2, 3
%     trellis = poly2trellis(3, [5 7]);
%     num_states = 4;
% 
%     % Frame Length
%     L_frame = 100; 
% 
%     errs = 0; total_bits = 0; blocks = 0;
% 
%     while errs < target_err && blocks < max_blocks
%         % 1. Generate Bits
%         data = randi([0 1], L_frame, 1);
% 
%         % 2. Encode (Convolutional Code)
%         % Output is [N, 2] matrix. Col 1 -> Ant 1, Col 2 -> Ant 2
%         enc_bits = convenc(data, trellis);
%         tx_syms_idx = reshape(enc_bits, 2, []).'; % [L_frame x 2]
% 
%         % 3. Map to BPSK
%         X_frame = 1 - 2 * tx_syms_idx; % [L_frame x 2]
%         X_frame = X_frame.';           % [2 x L_frame]
% 
%         % 4. Channel
%         if strcmp(fading_type, 'fast')
%             H_store = (randn(2,2,L_frame) + 1i*randn(2,2,L_frame))/sqrt(2);
%         else
%             H_blk = (randn(2,2)+1i*randn(2,2))/sqrt(2);
%             H_store = repmat(H_blk, [1, 1, L_frame]);
%         end
%         N_noise = (randn(2,L_frame) + 1i*randn(2,L_frame))/sqrt(2);
% 
%         Y = zeros(2, L_frame);
%         for t = 1:L_frame
%             Y(:,t) = H_store(:,:,t) * X_frame(:,t) * scale + N_noise(:,t);
%         end
% 
%         % 5. MIMO Viterbi Decoder (Soft Input)
%         % Initialize Path Metrics (Log Domain: Minimize Distance)
%         path_metrics = inf(num_states, 1);
%         path_metrics(1) = 0; % Start at state 0
% 
%         % Store survivor paths: [State, Time]
%         survivor_paths = zeros(num_states, L_frame);
% 
%         % Next States and Outputs from Trellis structure for speed
%         % trellis.nextStates: [States x Inputs] -> [4 x 2]
%         % trellis.outputs:    [States x Inputs] -> [4 x 2] (Decimal 0-3)
% 
%         for t = 1:L_frame
%             y_t = Y(:,t);
%             H_t = H_store(:,:,t);
%             new_path_metrics = inf(num_states, 1);
% 
%             % For each current state
%             for s = 1:num_states
%                 if path_metrics(s) == inf, continue; end
% 
%                 % Try both input bits 0 and 1
%                 for input_bit = 0:1
%                     % MATLAB trellis states are 0-indexed in logic, but indices are 1-based
%                     next_s = trellis.nextStates(s, input_bit+1) + 1; % 1-based index
%                     out_dec = trellis.outputs(s, input_bit+1);       % Decimal output (0-3)
% 
%                     % Convert decimal output to 2 BPSK symbols
%                     % out_dec = 2*b1 + b2. b1->Ant1, b2->Ant2
%                     % Note: convenc output order depends on poly order [5 7]. 
%                     % [5] is MSB usually. Let's decode manually:
%                     b1 = floor(out_dec/2);
%                     b2 = mod(out_dec, 2);
% 
%                     x_expected = [1 - 2*b1; 1 - 2*b2]; % BPSK Mapping
% 
%                     % Calculate Branch Metric (Euclidean Distance squared)
%                     % BM = || y - H*x ||^2
%                     bm = norm(y_t - H_t * x_expected * scale)^2;
% 
%                     % ACS (Add-Compare-Select)
%                     if path_metrics(s) + bm < new_path_metrics(next_s)
%                         new_path_metrics(next_s) = path_metrics(s) + bm;
%                         survivor_paths(next_s, t) = s; % Store previous state (1-based)
%                     end
%                 end
%             end
%             path_metrics = new_path_metrics;
%         end
% 
%         % 6. Traceback
%         [~, best_final_state] = min(path_metrics);
%         curr_s = best_final_state;
%         decoded_bits = zeros(L_frame, 1);
% 
%         for t = L_frame:-1:1
%             prev_s = survivor_paths(curr_s, t);
% 
%             % Find which input bit caused transition prev_s -> curr_s
%             % trellis.nextStates(prev_s, 1) is input 0
%             % trellis.nextStates(prev_s, 2) is input 1
%             if (trellis.nextStates(prev_s, 1) + 1) == curr_s
%                 decoded_bits(t) = 0;
%             else
%                 decoded_bits(t) = 1;
%             end
%             curr_s = prev_s;
%         end
% 
%         errs = errs + sum(abs(data - decoded_bits));
%         total_bits = total_bits + L_frame;
%         blocks = blocks + 1;
%     end
%     ber = errs / total_bits;
% end