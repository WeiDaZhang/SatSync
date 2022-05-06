% Baseband Processing
function BasebandProcessing(IQBB, ChipLength, ChipRate, nOverSample, LFSR_POL)

N = length(IQBB);
t = 1/ChipRate/nOverSample/2 : 1/ChipRate/nOverSample : (N-1/2)/ChipRate/nOverSample;

%Doppler_Removal_seq = ones(size(t));
%IQBB_doppler_remove_seq = [0 IQBB(1 : ChipLength*nOverSample - 1)];
Costas_Product = zeros(size(t));
WTune_y_n = zeros(size(t));
Doppler_Removal = ones(size(t));

% Low Pass Filter
lpf_b = [0.028  0.053 0.071  0.053 0.028];
lpf_a = [1.000 -2.026 2.148 -1.159 0.279];


% 2nd Order Filter Parameters
% Costas_Loop_Gain = 1e7;
% loop_f1 = 5e5;
% loop_f2 = 1e6;
% coef_b0 = 1;
% coef_b1 = pi*loop_f2*(1/ChipRate/nOverSample) + loop_f2/loop_f1;
% coef_b2 = pi*loop_f2*(1/ChipRate/nOverSample) - loop_f2/loop_f1;
% coef_a1 = 1;
% coef_a2 = 0;

% 3rd Order Filter Parameters locks
% Costas_Loop_Gain = 1e7;
% loop_tau1 = 1e-4;
% loop_tau2 = 0.2e-5;

Costas_Loop_Gain = 1e7;
loop_tau1 = 1e-5;
loop_tau2 = 2e-6;
coef_b0 = 1;
coef_b1 = 2*(1/ChipRate/nOverSample)/(1/ChipRate/nOverSample + 2*loop_tau2);
coef_b2 = (1/ChipRate/nOverSample - 2*loop_tau2)/(1/ChipRate/nOverSample + 2*loop_tau2);
coef_a1 = 4*loop_tau1/(1/ChipRate/nOverSample + 2*loop_tau1);
coef_a2 = (1/ChipRate/nOverSample - 2*loop_tau1)/(1/ChipRate/nOverSample + 2*loop_tau1);


LFSR_INIT = [0 0 0 0 0 1];
LocalCodeReplica = LFSR_t(t, ChipRate, LFSR_INIT, LFSR_POL, ChipLength);
idx_CodeReplica = 1;
Integral_Early_seq = zeros(size(IQBB));
Integral_Prompt_seq = zeros(size(IQBB));
Integral_Late_seq = zeros(size(IQBB));
Integral_dot = zeros(size(IQBB));

CodeLockThreshold = 8e-3;

h1 = figure;
title('Code Sequence')
h2 = figure;
title('Doppler Signal')
h3 = figure;
title('PDFx_n')
h4 = figure;
title('Integral Dot Product')

for idx_moving = 1 : (N - ChipLength*nOverSample)
    % Latest sequence of ChipRate*nOverSample samples of IQBB
    IQBB_seq = IQBB(idx_moving : ChipLength*nOverSample + idx_moving - 1);
    
    % IQ Mix with updated doppler removal NCO
    IQBB_doppler_remove_seq = IQBB_seq .* Doppler_Removal(idx_moving : ChipLength*nOverSample + idx_moving - 1);
    
    % Half Band Filter
    IQBB_doppler_remove_seq = filter(lpf_b,lpf_a,IQBB_doppler_remove_seq);

    % Phase-Frequency Detector (PFD) of Costas Loop
    PDFx_n = ...
                                 sum(  real(IQBB_doppler_remove_seq(end)) ... 
                                   .*  imag(IQBB_doppler_remove_seq(end)));

    % IIR Loop Filter
    WTune_y_n(ChipLength*nOverSample + idx_moving) = (coef_a1 * WTune_y_n(ChipLength*nOverSample + idx_moving - 1) + ...
                                                      coef_a2 * WTune_y_n(ChipLength*nOverSample + idx_moving - 2) + ...
                                    Costas_Loop_Gain*(coef_b0 * PDFx_n + ...
                                                      coef_b1 * Costas_Product(ChipLength*nOverSample + idx_moving - 1) + ...
                                                      coef_b2 * Costas_Product(ChipLength*nOverSample + idx_moving - 2) ...
                                                  ));

    % PFD Output
    Costas_Product(ChipLength*nOverSample + idx_moving) = PDFx_n;

    % Doppler Removal NCO
    Doppler_Removal(ChipLength*nOverSample + idx_moving) = exp(1j*WTune_y_n(ChipLength*nOverSample + idx_moving)*t(ChipLength*nOverSample + idx_moving));

    % Baseband signal Correlate with Early/Prompt/Late Code
    Correlated_Early_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica : ChipLength*nOverSample + idx_CodeReplica - 1);
    Correlated_Prompt_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica + 1 : ChipLength*nOverSample + idx_CodeReplica);
    Correlated_Late_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica + 2 : ChipLength*nOverSample + idx_CodeReplica + 1);
    
    % Integrate Correlation Result
    Integral_Early_seq(idx_moving) = sum(Correlated_Early_seq);
    Integral_Prompt_seq(idx_moving) = sum(Correlated_Prompt_seq);
    Integral_Late_seq(idx_moving) = sum(Correlated_Late_seq);
    Integral_dot(idx_moving) = ...
        (real(Integral_Early_seq(idx_moving)) - real(Integral_Late_seq(idx_moving))) ...
                                           .*real(Integral_Prompt_seq(idx_moving)) + ...
        (imag(Integral_Early_seq(idx_moving)) - imag(Integral_Late_seq(idx_moving))) ...
                                           .*imag(Integral_Prompt_seq(idx_moving));

    % Chip Lock Detection
    if((Integral_dot(idx_moving) > CodeLockThreshold) && (idx_moving > 1))
        if(Integral_dot(idx_moving - 1) < CodeLockThreshold)
            
            plot_graphs(h1, h2, h3, h4, IQBB_doppler_remove_seq, LocalCodeReplica, ...
                     idx_CodeReplica, ChipLength, nOverSample, Doppler_Removal, ...
                     Costas_Product, WTune_y_n, Integral_dot)
            disp('Peak-Start')
        end
        idx_CodeReplica = idx_CodeReplica + 1;
    end
    if((Integral_dot(idx_moving) < CodeLockThreshold) && (idx_moving > 1))
        if(Integral_dot(idx_moving - 1) > CodeLockThreshold)

            plot_graphs(h1, h2, h3, h4, IQBB_doppler_remove_seq, LocalCodeReplica, ...
                     idx_CodeReplica, ChipLength, nOverSample, Doppler_Removal, ...
                     Costas_Product, WTune_y_n, Integral_dot)
            disp('Peak-End')
        end
    end
    if(idx_moving == 0)
        plot_graphs(h1, h2, h3, h4, IQBB_doppler_remove_seq, LocalCodeReplica, ...
                     idx_CodeReplica, ChipLength, nOverSample, Doppler_Removal, ...
                     Costas_Product, WTune_y_n, Integral_dot)
        disp('Pause')
    end
end

plot_graphs(h1, h2, h3, h4, IQBB_doppler_remove_seq, LocalCodeReplica, ...
                     idx_CodeReplica, ChipLength, nOverSample, Doppler_Removal, ...
                     Costas_Product, WTune_y_n, Integral_dot)
disp('End')
end

function plot_graphs(h1, h2, h3, h4, IQBB_doppler_remove_seq, LocalCodeReplica, ...
                     idx_CodeReplica, ChipLength, nOverSample, Doppler_Removal, ...
                     Costas_Product, WTune_y_n, Integral_dot)
    figure(h1);
    plot(imag(IQBB_doppler_remove_seq));
    yyaxis right;
    plot(LocalCodeReplica(idx_CodeReplica + 1 : ChipLength*nOverSample + idx_CodeReplica));
    ylim([-2 2]);
    yyaxis left;
    title('Code Sequence');

    figure(h2);
    plot(real(Doppler_Removal));
    hold on;
    plot(imag(Doppler_Removal));
    title('Doppler Signal')

    figure(h3);
    plot(Costas_Product);
    yyaxis right;
    plot(WTune_y_n);
    yyaxis left;
    title('PDFx_n');

    figure(h4);
    plot((Integral_dot));
    title('Integral Dot Product');
end
