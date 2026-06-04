clear;
clc;
close all;

% System Parameters
alpha_pl = 3.0;
d_AR     = 1.0;
d_RB     = 1.0;
d_AE     = 0.8;
d_RE     = 0.8;
Omega_AR = d_AR^(-alpha_pl);
Omega_RB = d_RB^(-alpha_pl);
Omega_AE = d_AE^(-alpha_pl);
Omega_RE = d_RE^(-alpha_pl);

% Fixed SNR
SNR_dB_fixed = 10;
SNR_lin_fixed = 10^(SNR_dB_fixed/10);

% Secrecy Rate 
Rs_vec = 0:0.05:3.0;   

% Numerical Integration Settings
y_max   = 1e4;
rel_tol = 1e-8;
abs_tol = 1e-10;

% M values
M_scheme1 = [2, 3, 5];
M_scheme2 = [2, 3, 5];

% Plot Setup
figure('Color', 'w', 'Position', [100 100 800 550]);
hold on;

markers_s1 = {'-s', '-o', '-^'};
markers_s2 = {'--s', '--o', '--^'};
colors_s1  = lines(length(M_scheme1));
colors_s2  = lines(length(M_scheme2));

legendentries = {};
allEST = [];


% SCHEME 1: ZFB

fprintf('=== Scheme 1: ZFB ===\n');
for m_idx = 1:length(M_scheme1)
    M   = M_scheme1(m_idx);
    EST = zeros(size(Rs_vec));
    fprintf('Computing ZFB EST for M = %d...\n', M);

    for r_idx = 1:length(Rs_vec)
        Rs = Rs_vec(r_idx);
        A  = 2^(Rs) - 1;
        B  = 2^(Rs);

        g_RB = SNR_lin_fixed;
        g_AR = g_RB * (Omega_AR / Omega_RB);
        g_AE = g_RB * (Omega_AE / Omega_RB);
        g_RE = g_RB * (Omega_RE / Omega_RB);

        lam_AE = 1 / g_AE;
        lam_AR = 1 / g_AR;
        lam_RB = 1 / g_RB;
        mu     = g_RE / g_AE;

        % Eavesdropper PDF - ZFB 
        f_E = @(y) (lam_AE*(M-1)*mu .* exp(-lam_AE*y) ./ (1 + mu*y).^M) + ...
                   (lam_AE           .* exp(-lam_AE*y) ./ (1 + mu*y).^(M-1));

        % Receiver Outage Term S(y)
        S = @(y) exp(-lam_AR*(A + B*y)) .* arrayfun(@(yy) ...
            sum( (lam_AR*(A+B*yy)).^(0:M-1) ./ factorial(0:M-1) ), y);

        % T1
        T1 = integral(f_E, 0, y_max, 'RelTol', rel_tol, 'AbsTol', abs_tol);

        % T2
        T2 = integral(@(y) -S(y) .* f_E(y), 0, y_max, ...
                      'RelTol', rel_tol, 'AbsTol', abs_tol);

        % T3
        T3 = 0;
        for j = 1:M
            coeff_j = nchoosek(M,j) * (-1)^j * exp(-j*lam_RB*A);
            I3j     = integral(@(y) exp(-j*B*lam_RB*y) .* f_E(y), 0, y_max, ...
                               'RelTol', rel_tol, 'AbsTol', abs_tol);
            T3      = T3 + coeff_j * I3j;
        end

        % T4
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                Psi_jk   = (j*lam_RB + lam_AR);
                coeff_jk = nchoosek(M,j) * (-1)^j * (lam_AR^k / factorial(k)) * ...
                           exp(-Psi_jk * A);
                I4jk     = integral(@(y) (A + B*y).^k .* exp(-Psi_jk*B*y) .* ...
                           f_E(y), 0, y_max, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                T4       = T4 - coeff_jk * I4jk;
            end
        end

        SOP = T1 + T2 + T3 + T4;
        SOP = max(0, min(1, SOP));   % clamp to [0,1]

        % Effective Secrecy Throughput
        EST(r_idx) = Rs * (1 - SOP);
    end

    allEST = [allEST, EST(:)];
    plot(Rs_vec, EST, markers_s1{m_idx}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors_s1(m_idx,:), 'Color', colors_s1(m_idx,:), ...
        'MarkerIndices', 1:4:length(Rs_vec));
    legendentries{end+1} = sprintf('ZFB, M = %d', M);
end


% NAN SCHEME

fprintf('\n=== NAN Scheme ===\n');
for m_idx = 1:length(M_scheme2)
    M   = M_scheme2(m_idx);
    EST = zeros(size(Rs_vec));
    fprintf('Computing NAN EST for M = %d...\n', M);

    for r_idx = 1:length(Rs_vec)
        Rs = Rs_vec(r_idx);
        A  = 2^(Rs) - 1;
        B  = 2^(Rs);

        rho  = SNR_lin_fixed;
        g_AR = rho * Omega_AR;
        g_RB = rho * Omega_RB;
        g_AE = rho * Omega_AE;
        g_RE = rho * Omega_RE;

        beta   = g_RE / ((M-1) * g_AE);
        alpha  = g_RE / g_AE;
        lam_AE = 1 / g_AE;
        lam_AR = 1 / g_AR;
        lam_RB = 1 / g_RB;

        % Eavesdropper PDF - NAN 
        f_E = @(y) (lam_AE ./ (1 + beta.*y).^(M-1) + ...
                    alpha .* lam_AE .* y ./ (1 + beta.*y).^M) .* exp(-lam_AE.*y);

        % T1
        T1 = integral(f_E, 0, Inf, 'RelTol', rel_tol, 'AbsTol', abs_tol);

        % T2
        T2 = 0;
        for k = 0:M-1
            integrand_T2 = @(y) ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                                exp(-lam_AR*(A + B.*y)) .* f_E(y);
            T2 = T2 - integral(integrand_T2, 0, Inf, 'RelTol', rel_tol, 'AbsTol', abs_tol);
        end

        % T3
        T3 = 0;
        for j = 1:M
            coeff_j      = nchoosek(M,j) * (-1)^j;
            integrand_T3 = @(y) coeff_j .* exp(-j * lam_RB * (A + B.*y)) .* f_E(y);
            T3           = T3 + integral(integrand_T3, 0, Inf, 'RelTol', rel_tol, 'AbsTol', abs_tol);
        end

        % T4
        T4 = 0;
        for j = 1:M
            for k = 0:M-1
                coeff_jk     = nchoosek(M,j) * (-1)^j;
                integrand_T4 = @(y) coeff_jk .* exp(-j * lam_RB * (A + B.*y)) .* ...
                                    ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                                    exp(-lam_AR*(A + B.*y)) .* f_E(y);
                T4 = T4 - integral(integrand_T4, 0, Inf, 'RelTol', rel_tol, 'AbsTol', abs_tol);
            end
        end

        SOP = T1 + T2 + T3 + T4;
        SOP = max(0, min(1, SOP));   % clamp to [0,1]

        % Effective Secrecy Throughput
        EST(r_idx) = Rs * (1 - SOP);
    end

    allEST = [allEST, EST(:)];
    plot(Rs_vec, EST, markers_s2{m_idx}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors_s2(m_idx,:), 'Color', colors_s2(m_idx,:), ...
        'MarkerIndices', 1:4:length(Rs_vec));
    legendentries{end+1} = sprintf('NAN, M = %d', M);
end


grid on; grid minor;
xlabel('Secrecy Rate R_s (bits/s/Hz)', 'FontSize', 12);
ylabel('Effective Secrecy Throughput (EST)', 'FontSize', 12);
title(sprintf('EST vs Secrecy Rate (Fixed SNR = %d dB)', SNR_dB_fixed), 'FontSize', 13);
legend(legendentries, 'Location', 'NorthEast', 'FontSize', 10, 'NumColumns', 2);
xlim([Rs_vec(1), Rs_vec(end)]);
ylim([0, max(allEST(:)) * 1.1]);

hold off;
fprintf('\nPlotting Complete.\n');

