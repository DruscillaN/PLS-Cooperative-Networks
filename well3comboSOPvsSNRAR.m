clear;
clc;
close all;

% System Parameters
Rs       = 0.5;
A        = 2^(Rs) - 1;
B        = 2^(Rs);
alpha_pl = 3.0;
d_AR     = 1.0;
d_RB     = 1.0;
d_AE     = 0.8;
d_RE     = 0.8;
Omega_AR = d_AR^(-alpha_pl);
Omega_RB = d_RB^(-alpha_pl);
Omega_AE = d_AE^(-alpha_pl);
Omega_RE = d_RE^(-alpha_pl);

% SNR_AR range (x-axis)
SNR_AR_dB  = 0:2:40;
SNR_AR_lin = 10.^(SNR_AR_dB/10);

% Numerical Integration Settings
y_max   = 1e4;
rel_tol = 1e-8;
abs_tol = 1e-10;

% M values per scheme
M_scheme1 = [2, 3, 5];
M_scheme2 = [2, 3, 5];

% Plot Setup
figure('Color', 'w', 'Position', [100 100 800 550]);
hold on;

markers_s1 = {'-s', '-o', '-^', '-d'};
markers_s2 = {'--s', '--o', '--^'};
colors_s1  = lines(length(M_scheme1));
colors_s2  = lines(length(M_scheme2));

legendentries = {};
allSOP = [];


% SCHEME 1: ZFB

fprintf('=== Scheme 1 (ZFB) ===\n');
for m_idx = 1:length(M_scheme1)
    M   = M_scheme1(m_idx);
    SOP = zeros(size(SNR_AR_lin));
    fprintf('Computing Scheme 1 SOP for M = %d...\n', M);

    for idx = 1:length(SNR_AR_lin)
        % Derive all link SNRs from SNR_AR using path loss ratios
        g_AR = SNR_AR_lin(idx);                      % This IS SNR_AR
        g_RB = g_AR * (Omega_RB / Omega_AR);
        g_AE = g_AR * (Omega_AE / Omega_AR);
        g_RE = g_AR * (Omega_RE / Omega_AR);

        lam_AE = 1 / g_AE;
        lam_AR = 1 / g_AR;
        lam_RB = 1 / g_RB;
        mu     = g_RE / g_AE;

        % Eavesdropper PDF - Scheme 1 ZFB (UNCHANGED)
        f_E = @(y) (lam_AE*(M-1)*mu .* exp(-lam_AE*y) ./ (1 + mu*y).^M) + ...
                   (lam_AE           .* exp(-lam_AE*y) ./ (1 + mu*y).^(M-1));

        % Receiver Outage Term S(y) (UNCHANGED)
        S = @(y) exp(-lam_AR*(A + B*y)) .* arrayfun(@(yy) ...
            sum( (lam_AR*(A+B*yy)).^(0:M-1) ./ factorial(0:M-1) ), y);

        % T1 (UNCHANGED)
        T1 = integral(f_E, 0, y_max);

        % T2 (UNCHANGED)
        T2 = integral(@(y) -S(y) .* f_E(y), 0, y_max);

        % T3 (UNCHANGED)
        T3 = 0;
        for j = 1:M
            coeff_j = nchoosek(M,j) * (-1)^j * exp(-j*lam_RB*A);
            I3j     = integral(@(y) exp(-j*B*lam_RB*y) .* f_E(y), 0, y_max);
            T3      = T3 + coeff_j * I3j;
        end

        % T4 (UNCHANGED)
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                Psi_jk   = (j*lam_RB + lam_AR);
                coeff_jk = nchoosek(M,j) * (-1)^j * (lam_AR^k / factorial(k)) * ...
                           exp(-Psi_jk * A);
                I4jk     = integral(@(y) (A + B*y).^k .* exp(-Psi_jk*B*y) .* ...
                           f_E(y), 0, y_max);
                T4       = T4 - coeff_jk * I4jk;
            end
        end

        SOP(idx) = T1 + T2 + T3 + T4;
    end

    allSOP = [allSOP, SOP];
    semilogy(SNR_AR_dB, SOP, markers_s1{m_idx}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors_s1(m_idx,:), 'Color', colors_s1(m_idx,:));
    legendentries{end+1} = sprintf('Scheme 1 (ZFB), M = %d', M);
end


% SCHEME 2: NAN

fprintf('\n=== NAN Scheme ===\n');
for m_idx = 1:length(M_scheme2)
    M   = M_scheme2(m_idx);
    SOP = zeros(size(SNR_AR_lin));
    fprintf('Computing NAN SOP for M = %d...\n', M);

    for idx = 1:length(SNR_AR_lin)
        % Derive all link SNRs from SNR_AR using path loss ratios
        g_AR = SNR_AR_lin(idx);                      % This IS SNR_AR
        g_RB = g_AR * (Omega_RB / Omega_AR);
        g_AE = g_AR * (Omega_AE / Omega_AR);
        g_RE = g_AR * (Omega_RE / Omega_AR);

        beta   = g_RE / ((M-1) * g_AE);
        alpha  = g_RE / g_AE;
        lam_AE = 1/g_AE;
        lam_AR = 1/g_AR;
        lam_RB = 1/g_RB;

        % Eavesdropper PDF - NAN Scheme 
        f_E = @(y) (lam_AE ./ (1 + beta.*y).^(M-1) + ...
                    alpha .* lam_AE .* y ./ (1 + beta.*y).^M) .* exp(-lam_AE.*y);

        % T1 
        T1 = integral(f_E, 0, Inf);

        % T2 
        T2 = 0;
        for k = 0:M-1
            integrand_T2 = @(y) ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                                exp(-lam_AR*(A + B.*y)) .* f_E(y);
            T2 = T2 - integral(integrand_T2, 0, Inf);
        end

        % T3 
        T3 = 0;
        for j = 1:M
            coeff_j      = nchoosek(M,j) * (-1)^j;
            integrand_T3 = @(y) coeff_j .* exp(-j * lam_RB * (A + B.*y)) .* f_E(y);
            T3           = T3 + integral(integrand_T3, 0, Inf);
        end

        % T4 
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                coeff_jk     = nchoosek(M,j) * (-1)^j;
                integrand_T4 = @(y) coeff_jk .* exp(-j * lam_RB * (A + B.*y)) .* ...
                                    ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                                    exp(-lam_AR*(A + B.*y)) .* f_E(y);
                T4 = T4 - integral(integrand_T4, 0, Inf);
            end
        end

        SOP(idx) = T1 + T2 + T3 + T4;
        fprintf('SNR_AR=%-6.1f dB | T1=%-10.4f T2=%-10.4f T3=%-10.4f T4=%-10.4f SOP=%-12.4e\n', ...
                SNR_AR_dB(idx), T1, T2, T3, T4, SOP(idx));
    end

    allSOP = [allSOP, SOP];
    semilogy(SNR_AR_dB, SOP, markers_s2{m_idx}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors_s2(m_idx,:), 'Color', colors_s2(m_idx,:));
    legendentries{end+1} = sprintf('NAN, M = %d', M);
end

grid on; grid minor;
xlabel('Average SNR_{AR} (dB)', 'FontSize', 12);
ylabel('Secrecy Outage Probability (SOP)', 'FontSize', 12);
title(['SOP vs SNR_{AR}  (R_s = ', num2str(Rs), ' bits/s/Hz)'], 'FontSize', 13);
legend(legendentries, 'Location', 'SouthWest', 'FontSize', 10, 'NumColumns', 2);

ax = gca;
set(ax, 'YScale', 'log', 'YMinorGrid', 'on');
ymin_exp  = floor(log10(min(allSOP(allSOP > 0))));
ymax_exp  = ceil(log10(max(allSOP(allSOP > 0))));
ytick_exp = ymin_exp:ymax_exp;
yticks(10.^ytick_exp);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), ytick_exp, 'UniformOutput', false));
ylim([10^ymin_exp, 10^ymax_exp]);

hold off;
fprintf('\nPlotting Complete.\n');
