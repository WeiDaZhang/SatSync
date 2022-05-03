% Baseband Processing
function BasebandProcessing(IQBB, ChipLength, ChipRate, nOverSample, LFSR_POL)

N = length(IQBB);
t = 1/ChipRate/nOverSample/2 : 1/ChipRate/nOverSample : (N-1/2)/ChipRate/nOverSample;

%Doppler_Removal_seq = ones(size(t));
IQBB_doppler_remove_seq = [0 IQBB(1 : ChipLength*nOverSample - 1)];
Costas_Product = zeros(size(t));

% 2nd Order Filter Parameters
% Costas_Loop_Gain = 1e9;
% loop_f1 = 5e5;
% loop_f2 = 1e5;
% coef1 = pi*loop_f2*(1/ChipRate/nOverSample) + loop_f2/loop_f1;
% coef2 = pi*loop_f2*(1/ChipRate/nOverSample) - loop_f2/loop_f1;

% 3rd Order Filter Parameters
Costas_Loop_Gain = 1e9;
loop_tau1 = 1e-3;
loop_tau2 = 0.2e-4;
coef_b0 = 1;
coef_b1 = 2*(1/ChipRate/nOverSample)/(1/ChipRate/nOverSample + 2*loop_tau2);
coef_b2 = (1/ChipRate/nOverSample - 2*loop_tau2)/(1/ChipRate/nOverSample + 2*loop_tau2);
coef_a1 = 4*loop_tau1/(1/ChipRate/nOverSample + 2*loop_tau1);
coef_a2 = (1/ChipRate/nOverSample - 2*loop_tau1)/(1/ChipRate/nOverSample + 2*loop_tau1);

WTune_y_n = zeros(size(t));

LFSR_INIT = [0 0 0 0 0 1];
LocalCodeReplica = LFSR_t(t, ChipRate, LFSR_INIT, LFSR_POL, ChipLength);
idx_CodeReplica = 1;
Integral_Early_seq = zeros(size(IQBB));
Integral_Prompt_seq = zeros(size(IQBB));
Integral_Late_seq = zeros(size(IQBB));
Integral_dot = zeros(size(IQBB));
for idx_moving = 1 : (N - ChipLength*nOverSample)
    % Latest sequence of ChipRate*nOverSample samples of IQBB
    IQBB_seq = IQBB(idx_moving : ChipLength*nOverSample + idx_moving - 1);

    PDFx_n = ...
                                 sum( [ real(IQBB_doppler_remove_seq(2 : ...
                                                                     ChipLength*nOverSample) ...
                                             ) ... 
                                        real(IQBB_seq(end)) ...
                                      ] ... 
                                   .* [ imag(IQBB_doppler_remove_seq(2 : ...
                                                                     ChipLength*nOverSample) ...
                                             ) ... 
                                        imag(IQBB_seq(end)) ...
                                      ] ...
                                    );
    WTune_y_n(ChipLength*nOverSample + idx_moving) = (coef_a1 * WTune_y_n(ChipLength*nOverSample + idx_moving - 1) + ...
                                                      coef_a2 * WTune_y_n(ChipLength*nOverSample + idx_moving - 2) + ...
                                    Costas_Loop_Gain*(coef_b0 * PDFx_n + ...
                                                      coef_b1 * Costas_Product(ChipLength*nOverSample + idx_moving - 1) + ...
                                                      coef_b2 * Costas_Product(ChipLength*nOverSample + idx_moving - 2) ...
                                                  ));

    Costas_Product(ChipLength*nOverSample + idx_moving) = PDFx_n;

%    Doppler_Removal_seq(ChipLength*nOverSample - 1 + idx_moving) = exp(1j*WTune_y_n(ChipLength*nOverSample + idx_moving)*t(ChipLength*nOverSample - 1 + idx_moving));
    Doppler_Removal_seq = exp(1j*WTune_y_n(ChipLength*nOverSample + idx_moving)*t(idx_moving : ChipLength*nOverSample + idx_moving - 1));
%    IQBB_doppler_remove_seq = IQBB_seq .* Doppler_Removal_seq(idx_moving : ChipLength*nOverSample + idx_moving - 1);
    IQBB_doppler_remove_seq = IQBB_seq .* Doppler_Removal_seq;
    
    %% Filter for removing the image components of doppler removal is needed
    
    %%
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
    if(Integral_dot(idx_moving) > 2.5e-4)
        figure;
        plot(real(IQBB_doppler_remove_seq))
        yyaxis right;
        plot(LocalCodeReplica(idx_CodeReplica + 1 : ChipLength*nOverSample + idx_CodeReplica))
        yyaxis left;
        figure;plot(real(Doppler_Removal_seq))
        hold on;plot(imag(Doppler_Removal_seq))
        disp('Peak')
        idx_CodeReplica = idx_CodeReplica + 1;
    end
    if(idx_moving == 0)
        figure;
        plot(real(IQBB_doppler_remove_seq))
        yyaxis right;
        plot(LocalCodeReplica(idx_CodeReplica + 1 : ChipLength*nOverSample + idx_CodeReplica))
        yyaxis left;
        figure;plot(real(Doppler_Removal_seq))
        hold on;plot(imag(Doppler_Removal_seq))
        disp('Stopped')
    end
end
figure;plot(Costas_Product)
title('PDFx_n')
yyaxis right
plot(WTune_y_n)

figure;plot((Integral_dot));
title('Integral Dot Product')
