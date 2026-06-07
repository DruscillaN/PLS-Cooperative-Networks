clear; 
clc; 
close all;


M_array  = [2, 3, 5];
Rs_array = 0:0.2:4;   % Secrecy rate range (bps/Hz) — 8 points
SNR_dB   = 20;

alpha = 3;
d_AR  = 1.0;
d_RB  = 1.0;
d_AE  = 1.5;
d_RE  = 0.8;

gamma_SR = 10^(SNR_dB/10) / d_SR^alpha;
gamma_SE = 10^(10/10)     / d_SE^alpha;
gamma_RE = 10^(10/10)     / d_RE^alpha;

fprintf('System Parameters\n');
fprintf('SNR = %d dB | alpha = %.1f\n', SNR_dB, alpha);
fprintf('d_SR=%.1f | d_SE=%.1f | d_RE=%.1f\n', d_SR, d_SE, d_RE);
fprintf('gamma_SR=%.4f | gamma_SE=%.4f | gamma_RE=%.4f\n\n', gamma_SR, gamma_SE, gamma_RE);


ST_ZFB_all = zeros(length(M_array), length(Rs_array));
ST_NAN_all = zeros(length(M_array), length(Rs_array));


colors = lines(length(M_array));
figure('Color','w','Position',[100 100 820 600]);
hold on;
grid on;

fprintf('Computing Secrecy Throughput vs Rs...\n');


for mm = 1:length(M_array)
    M = M_array(mm);
    ST_ZFB = zeros(size(Rs_array));
    ST_NAN = zeros(size(Rs_array));

    for ri = 1:length(Rs_array)
        Rs  = Rs_array(ri);
        gth = 2^(2*Rs) - 1;
        B   = 2^(2*Rs);

        %% ZFB CLOSED-FORM
        gamma_E_eff_ZFB = gamma_SE + gamma_RE / M;
        lam = B/gamma_SR + 1/gamma_E_eff_ZFB;
        s = 0;
        for k = 0:M-1
            for n = 0:k
                s = s + nchoosek(k,n) * gth^(k-n) * B^n * factorial(n) / ...
                    (factorial(k) * gamma_SR^k * gamma_E_eff_ZFB * lam^(n+1));
            end
        end
        SOP_ZFB    = 1 - exp(-gth/gamma_SR) * s;
        ST_ZFB(ri) = (1 - SOP_ZFB) * Rs;

        %% NAN CLOSED-FORM
        gamma_E_eff_NAN = gamma_SE + gamma_RE;
        lam_n = B/gamma_SR + 1/gamma_E_eff_NAN;
        s = 0;
        for k = 0:M
            for n = 0:k
                s = s + nchoosek(k,n) * gth^(k-n) * B^n * factorial(n) / ...
                    (factorial(k) * gamma_SR^k * gamma_E_eff_NAN * lam_n^(n+1));
            end
        end
        SOP_NAN    = 1 - exp(-gth/gamma_SR) * s;
        ST_NAN(ri) = (1 - SOP_NAN) * Rs;
    end

  
    ST_ZFB_all(mm, :) = ST_ZFB;
    ST_NAN_all(mm, :) = ST_NAN;

    fprintf('M=%d | ST_ZFB max=%.4f | ST_NAN max=%.4f\n', M, max(ST_ZFB), max(ST_NAN));

 
    plot(Rs_array, ST_ZFB, '-o', 'Color', colors(mm,:), 'LineWidth', 2, ...
        'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(mm,:), ...
        'DisplayName', sprintf('ZFB (M=%d)', M));
    plot(Rs_array, ST_NAN, '--s', 'Color', colors(mm,:), 'LineWidth', 2, ...
        'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(mm,:), ...
        'DisplayName', sprintf('NAN (M=%d)', M));
end

%% Labels & Title
xlabel('Secrecy Rate Threshold R_s (bps/Hz)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Secrecy Throughput ST (bps/Hz)',       'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Secrecy Throughput vs R_s | \\gamma_{SR}=%d dB, \\alpha=%.1f', ...
      SNR_dB, alpha), 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10, 'NumColumns', 1);

%% Axis Formatting
ax = gca;
ax.XLim          = [min(Rs_array), max(Rs_array)];
ax.YLim          = [0, ceil(max(ST_ZFB_all(:))*10)/10 + 0.1];
ax.FontSize      = 12;
ax.GridAlpha     = 0.4;
ax.GridLineStyle = '--';
ax.MinorGridLineStyle = ':';
grid(ax, 'on');

%% Save
save('slot_ST.mat', ...
     'ST_ZFB_all', 'ST_NAN_all', ...
     'Rs_array', 'M_array', 'SNR_dB', 'alpha', ...
     'd_SR', 'd_SE', 'd_RE', ...
     'gamma_SR', 'gamma_SE', 'gamma_RE');

fprintf('\nSaved to slot_ST.mat\n');
fprintf('Matrix dimensions: ST_ZFB_all is [%d x %d]\n', size(ST_ZFB_all));
