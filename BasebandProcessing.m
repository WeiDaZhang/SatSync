% Baseband Processing
function BasebandProcessing(IQBB, ChipLength, ChipRate, nOverSample, LFSR_POL)

Doppler_Removal_seq = ones(1, ChipLength*nOverSample);
N = length(IQBB);
t = 0 : 1/ChipRate/nOverSample : (N-1)/ChipRate/nOverSample;
LFSR_INIT = [0 0 0 0 0 1];
LocalCodeReplica = LFSR_t(t, ChipRate, LFSR_INIT, LFSR_POL);
idx_CodeReplica = 1;
Integral_Early_seq = zeros(size(IQBB));
Integral_Prompt_seq = zeros(size(IQBB));
Integral_Late_seq = zeros(size(IQBB));
Integral_dot = zeros(size(IQBB));
for idx_moving = 1 : (N - ChipLength*nOverSample)
    % Latest sequence of ChipRate*nOverSample samples of IQBB
    IQBB_seq = IQBB(1 + idx_moving : ChipLength*nOverSample + idx_moving);
    IQBB_doppler_remove_seq = IQBB_seq .* Doppler_Removal_seq;
    Correlated_Early_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica : ChipLength*nOverSample + idx_CodeReplica - 1);
    Correlated_Prompt_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica + 1 : ChipLength*nOverSample + idx_CodeReplica);
    Correlated_Late_seq = IQBB_doppler_remove_seq .* LocalCodeReplica(idx_CodeReplica + 2 : ChipLength*nOverSample + idx_CodeReplica + 1);
    Integral_Early_seq(idx_moving) = sum(Correlated_Early_seq);
    Integral_Prompt_seq(idx_moving) = sum(Correlated_Prompt_seq);
    Integral_Late_seq(idx_moving) = sum(Correlated_Late_seq);
    Integral_dot(idx_moving) = ...
        (real(Integral_Early_seq(idx_moving)) - real(Integral_Late_seq(idx_moving))) ...
                                           .*real(Integral_Prompt_seq(idx_moving)) + ...
        (imag(Integral_Early_seq(idx_moving)) - imag(Integral_Late_seq(idx_moving))) ...
                                           .*imag(Integral_Prompt_seq(idx_moving));
    if(Integral_dot(idx_moving) > 2.5e-5)
        disp('Peak')
    end
end
figure;plot(abs(Integral_Early_seq));
hold on;
plot(abs(Integral_Prompt_seq));
plot(abs(Integral_Late_seq));

figure;plot(abs(Integral_dot));
