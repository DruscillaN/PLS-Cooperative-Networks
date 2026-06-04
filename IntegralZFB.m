clear; 
clc; 
close all;

% System Parameters
M       = 3;        % Number of relay antennas 
Rs      = 0.5;      % Target secrecy rate 

A = 2^(Rs) - 1;   
B = 2^(Rs);       

alpha_pl = 3;     
d_AR = 1.0; 
d_RB = 1.0; 
d_AE = 1.5; 
d_RE = 0.8;

Omega_AR = d_AR^(-alpha_pl);
Omega_RB = d_RB^(-alpha_pl);
Omega_AE = d_AE^(-alpha_pl);
Omega_RE = d_RE^(-alpha_pl);

% Sweep Range for SNR_RD (Relay-Destination Link)
SNR_RD_dB  = 0:2:40; 
SNR_RD_lin = 10.^(SNR_RD_dB/10);

% Numerical Integration Settings
y_max   = 1e4;
rel_tol = 1e-8;
abs_tol = 1e-10;

% Preallocate
SOP = zeros(size(SNR_RD_lin));

fprintf('Computing SOP vs SNR_RD for M = %d...\n', M);

%Main Loop over SNR_RD
for idx = 1:length(SNR_RD_lin)
    
    g_RB = SNR_RD_lin(idx);
    g_AR = g_RB * (Omega_AR / Omega_RB); 
    g_AE = g_RB * (Omega_AE / Omega_RB);
    g_RE = g_RB * (Omega_RE / Omega_RB);

    % Rate parameters
    lam_AE = 1 / g_AE;
    lam_AR = 1 / g_AR;
    lam_RB = 1 / g_RB;
    mu     = g_RE / g_AE;    

    % Eavesdropper PDF
    f_E = @(y) (lam_AE*(M-1)*mu .* exp(-lam_AE*y) ./ (1 + mu*y).^M) + ...
               (lam_AE           .* exp(-lam_AE*y) ./ (1 + mu*y).^(M-1));

    % Receiver Outage Term S(y) 
    S = @(y) exp(-lam_AR*(A + B*y)) .*arrayfun(@(yy)...
        sum( (lam_AR*(A+B*yy)).^(0:M-1) ./ factorial(0:M-1) ), y);

    % T1: PDF Normalization
    T1 = integral(f_E, 0, y_max, 'RelTol', rel_tol, 'AbsTol', abs_tol);

    % T2: Source-Relay Outage
    T2 = integral(@(y) -S(y) .* f_E(y), 0, y_max, 'RelTol', ...
        rel_tol, 'AbsTol', abs_tol);

    % T3: Destination Diversity
    T3 = 0;
    for j = 1:M
        coeff_j = nchoosek(M,j) * (-1)^j * exp(-j*lam_RB*A);
        I3j = integral(@(y) exp(-j*B*lam_RB*y) .* f_E(y), 0, y_max, ...
                       'RelTol', rel_tol, 'AbsTol', abs_tol);
        T3 = T3 + coeff_j * I3j;
    end

    % T4: Joint Outage (Cross-term)
    T4 = 0;
    for j = 1:M
        for k = 0:M-1
            Psi_jk = (j*lam_RB + lam_AR);
            coeff_jk = nchoosek(M,j) * (-1)^j * (lam_AR^k / factorial(k)) *...
                exp(-Psi_jk * A);
            I4jk = integral(@(y) (A + B*y).^k .* exp(-Psi_jk*B*y) .* ...
                f_E(y), 0, y_max,'RelTol', rel_tol, 'AbsTol', abs_tol);
            T4 = T4 - coeff_jk * I4jk;
        end
    end

    % Final SOP sum
    SOP(idx) = T1 + T2 + T3 + T4;
end

% Plotting Results
figure('Color', 'w', 'Position', [100 100 700 500]);
semilogy(SNR_RD_dB, SOP, 'b-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on; grid minor;
xlabel('Average Relay-Destination SNR \gamma_{RD} (dB)', 'FontSize', 12);
ylabel('Secrecy Outage Probability (SOP)', 'FontSize', 12);
title(['SOP vs \gamma_{RD} (M = ', num2str(M), ', R_s = ', num2str(Rs), ' bits/s/Hz)']);
legend(['Numerical Integration (M=', num2str(M), ')'], 'Location', 'SouthWest');

fprintf('Plotting Complete.\n');