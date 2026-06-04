clear; 
clc; 
close all;

%System Parameters
M       = 3;        % Number of antennas
Rs      = 0.5;      % Target secrecy rate 

A = 2^(Rs) - 1;   
B = 2^(Rs);       

% Path loss parameters
alpha_pl = 3.0;     
d_AR = 1.0; 
d_RB = 1.0; 
d_AE = 1.5; 
d_RE = 0.8;

Omega_AR = d_AR^(-alpha_pl); 
Omega_RB = d_RB^(-alpha_pl);
Omega_AE = d_AE^(-alpha_pl); 
Omega_RE = d_RE^(-alpha_pl);

% SNR range
SNR_dB  = 0:2:40; 
SNR_lin = 10.^(SNR_dB/10);

% Preallocate
SOP = zeros(size(SNR_lin));

fprintf('Starting NAN SOP Computation...\n');
fprintf('%-10s %-12s %-12s %-12s %-12s %-12s\n', 'SNR(dB)', 'T1', 'T2', 'T3', 'T4', 'SOP');

% Main Loop
for idx = 1:length(SNR_lin)
    rho = SNR_lin(idx);
    
    % Link average SNRs
    g_AR = rho * Omega_AR; 
    g_RB = rho * Omega_RB;
    g_AE = rho * Omega_AE; 
    g_RE = rho * Omega_RE;
    
    % Constants for NAN PDF
    beta   = g_RE / ((M-1) * g_AE); 
    alpha  = g_RE / g_AE;           
    lam_AE = 1/g_AE; lam_AR = 1/g_AR; lam_RB = 1/g_RB;

    % Define Eavesdropper PDF 
    f_E = @(y) (lam_AE ./ (1 + beta.*y).^(M-1) + ...
                alpha .* lam_AE .* y ./ (1 + beta.*y).^M) .* exp(-lam_AE.*y);

    % Term T1: Baseline
    T1 = integral(f_E, 0, Inf);

    % Term T2: Source-Relay Outage 
    T2 = 0;
    for k = 0:M-1
        % Integrand for specific k: (gamma^k / k!) * exp(-gamma/g_AR) * f_E(y)
        integrand_T2 = @(y) ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                            exp(-lam_AR*(A + B.*y)) .* f_E(y);
        T2 = T2 - integral(integrand_T2, 0, Inf);
    end

    % Term T3: Relay-Destination Diversity 
    T3 = 0;
    for j = 1:M
        % Integrand for specific j: C(M,j) * (-1)^j * exp(-j*gamma/g_RB) * f_E(y)
        coeff_j = nchoosek(M,j) * (-1)^j;
        integrand_T3 = @(y) coeff_j .* exp(-j * lam_RB * (A + B.*y)) .* f_E(y);
        T3 = T3 + integral(integrand_T3, 0, Inf);
    end

    % Term T4: Cross-Product 
    T4 = 0;
    for j = 1:M
        for k = 0:M-1
            % Combined Outage: C(M,j) * (-1)^j * exp(-j*gamma/g_RB) * S(y) * f_E(y)
            coeff_jk = nchoosek(M,j) * (-1)^j;
            integrand_T4 = @(y) coeff_jk .* exp(-j * lam_RB * (A + B.*y)) .* ...
                                ((lam_AR*(A + B.*y)).^k ./ factorial(k)) .* ...
                                exp(-lam_AR*(A + B.*y)) .* f_E(y);
            T4 = T4 - integral(integrand_T4, 0, Inf);
        end
    end

    % Result for current SNR
    SOP(idx) = T1 + T2 + T3 + T4;
    fprintf('%-10.1f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4e\n', ...
            SNR_dB(idx), T1, T2, T3, T4, SOP(idx));
end

% Plotting
figure('Color', 'w');
semilogy(SNR_dB, SOP, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on; grid minor;
xlabel('Average Link SNR (dB)');
ylabel('Secrecy Outage Probability (SOP)');
title(['NAN Scheme SOP (M = ', num2str(M), ', R_s = ', num2str(Rs), ')']);