clear;
clc;
close all;

%% 1. System Parameters
SNR_fixed_dB  = 10;
SNR_fixed_lin = 10^(SNR_fixed_dB/10);

alpha_pl = 3.0;
d_AR = 1.0;
d_RB = 1.0;
d_AE = 0.8;
d_RE = 0.8;
Omega_AR = d_AR^(-alpha_pl);
Omega_RB = d_RB^(-alpha_pl);
Omega_AE = d_AE^(-alpha_pl);
Omega_RE = d_RE^(-alpha_pl);

% Sweep Range for Secrecy Rate Rs (bits/s/Hz)
Rs_vec = 0.1:0.5:2.6;


% M values per scheme
M_list = [2, 3, 5];

% Plot Setup
figure('Color', 'w', 'Position', [100 100 850 850]);
hold on;
colors_s1     = lines(length(M_list));
colors_s2     = lines(length(M_list));
legendentries = {};
allSOP        = [];

% Link average SNRs (Fixed based on SNR_fixed_dB)
g_AR = SNR_fixed_lin * Omega_AR;
g_RB = SNR_fixed_lin * Omega_RB;
g_AE = SNR_fixed_lin * Omega_AE;
g_RE = SNR_fixed_lin * Omega_RE;

lam_AE = 1 / g_AE;
lam_AR = 1 / g_AR;
lam_RB = 1 / g_RB;
mu     = g_RE / g_AE;

% Pre-allocate SOP storage matrices
SOP_ZFB_all = zeros(length(M_list), length(Rs_vec));
SOP_NAN_all = zeros(length(M_list), length(Rs_vec));

%% 2. SCHEME 1: ZFB (Eq. 24) vs Rs
fprintf('=== Scheme 1 (ZFB) vs Rs ===\n');
for m_idx = 1:length(M_list)
    M   = M_list(m_idx);
    SOP = zeros(size(Rs_vec));
    fprintf('Computing ZFB SOP for M = %d...\n', M);

    for idx = 1:length(Rs_vec)
        Rs_curr = Rs_vec(idx);
        % Thresholds from Eq. 30 (using 2*Rs for half-duplex)
        A = 2^(Rs_curr) - 1;
        B = 2^(Rs_curr);

        % Eavesdropper PDF: ZFB (Equation 24)
        f_E = @(y) (lam_AE*(M-1)*mu .* exp(-lam_AE*y) ./ (1 + mu*y).^M) + ...
                   (lam_AE           .* exp(-lam_AE*y) ./ (1 + mu*y).^(M-1));

        % Receiver Outage Term S(y)
        S = @(y) exp(-lam_AR*(A + B*y)) .* arrayfun(@(yy) ...
            sum( (lam_AR*(A+B*yy)).^(0:M-1) ./ factorial(0:M-1) ), y);

        T1 = integral(f_E, 0, Inf);
        T2 = integral(@(y) -S(y) .* f_E(y), 0, Inf);
        T3 = 0;
        for j = 1:M
            coeff_j = nchoosek(M,j) * (-1)^j * exp(-j*lam_RB*A);
            T3 = T3 + coeff_j * integral(@(y) exp(-j*B*lam_RB*y) .* f_E(y), 0, Inf);
        end
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                Psi_jk   = (j*lam_RB + lam_AR);
                coeff_jk = nchoosek(M,j) * (-1)^j * (lam_AR^k / factorial(k)) * exp(-Psi_jk * A);
                T4 = T4 - coeff_jk * integral(@(y) (A + B*y).^k .* exp(-Psi_jk*B*y) .* f_E(y), 0, Inf);
            end
        end
        SOP(idx) = T1 + T2 + T3 + T4;
    end

    SOP_ZFB_all(m_idx, :) = SOP;
    allSOP = [allSOP, SOP]; %#ok<AGROW>
    semilogy(Rs_vec, SOP, '-o', 'LineWidth', 2, 'MarkerSize', 5, 'Color', colors_s1(m_idx,:));
    legendentries{end+1} = sprintf('ZFB, M = %d', M); %#ok<SAGROW>
end

% Extract named per-M ZFB results
SOP_ZFB_M2 = SOP_ZFB_all(1, :);
SOP_ZFB_M3 = SOP_ZFB_all(2, :);
SOP_ZFB_M5 = SOP_ZFB_all(3, :);

%% 3. SCHEME 2: NAN (Eq. 25) vs Rs
fprintf('\n=== Scheme 2 (NAN) vs Rs ===\n');
for m_idx = 1:length(M_list)
    M   = M_list(m_idx);
    SOP = zeros(size(Rs_vec));
    fprintf('Computing NAN SOP for M = %d...\n', M);

    for idx = 1:length(Rs_vec)
        Rs_curr = Rs_vec(idx);
        A = 2^(Rs_curr) - 1;
        B = 2^(Rs_curr);

        beta  = g_RE / ((M-1) * g_AE);
        alpha = g_RE / g_AE;

        % Eavesdropper PDF: NAN (Equation 25)
        f_E = @(y) (lam_AE ./ (1 + beta.*y).^(M-1) + ...
                    alpha .* lam_AE .* y ./ (1 + beta.*y).^M) .* exp(-lam_AE.*y);

        T1 = integral(f_E, 0, Inf);
        T2 = 0;
        for k = 0:M-1
            T2 = T2 - integral(@(y) ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                 exp(-lam_AR*(A + B.*y)) .* f_E(y), 0, Inf);
        end
        T3 = 0;
        for j = 1:M
            coeff_j = nchoosek(M,j) * (-1)^j;
            T3 = T3 + coeff_j * integral(@(y) exp(-j * lam_RB * (A + B.*y)) .* f_E(y), 0, Inf);
        end
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                coeff_jk = nchoosek(M,j) * (-1)^j;
                T4 = T4 - coeff_jk * integral(@(y) exp(-j * lam_RB * (A + B.*y)) .* ...
                     ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* exp(-lam_AR*(A + B.*y)) .* f_E(y), 0, Inf);
            end
        end
        SOP(idx) = T1 + T2 + T3 + T4;
    end

    SOP_NAN_all(m_idx, :) = SOP;
    allSOP = [allSOP, SOP]; %#ok<AGROW>
    semilogy(Rs_vec, SOP, '--s', 'LineWidth', 2, 'MarkerSize', 5, 'Color', colors_s2(m_idx,:));
    legendentries{end+1} = sprintf('NAN, M = %d', M); %#ok<SAGROW>
end

% Extract named per-M NAN results
SOP_NAN_M2 = SOP_NAN_all(1, :);
SOP_NAN_M3 = SOP_NAN_all(2, :);
SOP_NAN_M5 = SOP_NAN_all(3, :);

%% 4. Formatting
grid on; grid minor;
xlabel('Target Secrecy Rate R_s (bits/s/Hz)', 'FontSize', 12);
ylabel('Secrecy Outage Probability (SOP)', 'FontSize', 12);
title(['SOP vs Secrecy Rate R_s (SNR = ', num2str(SNR_fixed_dB), ' dB)']);
legend(legendentries, 'Location', 'SouthEast', 'FontSize', 9, 'NumColumns', 2);

ax = gca;
ax.YScale = 'log';
ax.YTick  = 10.^(-12:1:0);
ytickformat('10^{%.0f}');
ylim([1e-7 1.1]);

hold off;
fprintf('\nPlotting Complete. X-axis represents Rs.\n');

%% 5. Save all variables to slot30.mat
save('slot30.mat', ...
    'SNR_fixed_dB', 'SNR_fixed_lin', ...
    'alpha_pl', ...
    'd_AR', 'd_RB', 'd_AE', 'd_RE', ...
    'Omega_AR', 'Omega_RB', 'Omega_AE', 'Omega_RE', ...
    'Rs_vec', ...
    'M_list', ...
    'g_AR', 'g_RB', 'g_AE', 'g_RE', ...
    'lam_AE', 'lam_AR', 'lam_RB', 'mu', ...
    'SOP_ZFB_M2', 'SOP_ZFB_M3', 'SOP_ZFB_M5', 'SOP_ZFB_all', ...
    'SOP_NAN_M2', 'SOP_NAN_M3', 'SOP_NAN_M5', 'SOP_NAN_all');

fprintf('\n=== All variables saved to slot30.mat ===\n');
fprintf('\nVariable Summary:\n');
fprintf('  Scalars   : SNR_fixed_dB, SNR_fixed_lin, alpha_pl, rel_tol, abs_tol\n');
fprintf('  Distances : d_AR, d_RB, d_AE, d_RE\n');
fprintf('  Gains     : Omega_AR, Omega_RB, Omega_AE, Omega_RE\n');
fprintf('  Link SNRs : g_AR, g_RB, g_AE, g_RE\n');
fprintf('  Exp rates : lam_AE, lam_AR, lam_RB, mu\n');
fprintf('  Rs axis   : Rs_vec [1x%d]\n', length(Rs_vec));
fprintf('  M set     : M_list\n');
fprintf('  ZFB SOP   : SOP_ZFB_M2, SOP_ZFB_M3, SOP_ZFB_M5, SOP_ZFB_all [%dx%d]\n', ...
        size(SOP_ZFB_all,1), size(SOP_ZFB_all,2));
fprintf('  NAN SOP   : SOP_NAN_M2, SOP_NAN_M3, SOP_NAN_M5, SOP_NAN_all [%dx%d]\n', ...
        size(SOP_NAN_all,1), size(SOP_NAN_all,2));