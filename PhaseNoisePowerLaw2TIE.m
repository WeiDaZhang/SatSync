% TimeIntervalError from Phase Noise Power Law
function delta_t = PhaseNoisePowerLaw2TIE(N, nominal_interval, S_phi_f0)
% N is the length of the time series

% nominal_interval (s) is the ideal interval

% S_phi_f0 is a vector of S_phi (phase noise dBc/Hz) for each f^beta slop
% intercepting f = 1Hz.

% beta = [-4,0] are random-walk-FM, flicker-FM, white-FM, flicker-PM, white-PM.

% size(s_phi_f0) = [5,1].

close all

if(mod(N,2))
    N = N + 1;
end
beta = -4:0;
S_phi = zeros(1,N);     % Phase (Amplitude) Spectrum Density (rad/sqrt(Hz))
bin_bandwidth = 1/nominal_interval / N; % delta frequency
for index = 1:5
    % S_phi input in rad^2/Hz is first convert to rad^2/(fft-bin)
    % ASD is then obtained by square root
    ASD = sqrt(S_phi_f0(index)/bin_bandwidth);

    % Slope of PSD is calculated with 5 beta values
    Slope = (1:N/2 - 1).^(beta(index)/2);

    % 1Hz interception point is fft-bin bandwidth to the power of sqrt(beta)
    Intercept = bin_bandwidth^(beta(index)/2);

    mag = ASD*Slope*Intercept.*randn(1,N/2 - 1);
    pha = 2*pi*rand(1,N/2 - 1);

    S_phi = S_phi + [0 mag.*exp(1j*pha) 0 conj(flip(mag.*exp(1j*pha)))];

    % For 2N point FFT, which is indexing [1:2N],
    % the first bin is the power middled in DC,
    % considering [0 : 2N-1]
    % the first bin contains (2N-0.5 to 2N][DC to 0.5)
    % For indexing [1:2N] the middle point is N+1
end

f = 0 : 1/((N -1) * nominal_interval) : 1/nominal_interval;
figure;
loglog(f, abs(S_phi.^2)*bin_bandwidth); 
grid on;
xlabel('Frequency / Hz');
ylabel('Phase Noise / rad^2/Hz')
title('Clock Difference in Phase Noise Spectrum Density')


phi = ((sqrt(N)*ifft(S_phi)*sqrt(bin_bandwidth)));

% For only white and flicker phase modulation
% Small difference should be found between sqrt(S_phi) and rms(phi)
% disp(sqrt(sum(S_phi_f0(4:5))))
% disp(rms(phi))
% disp(sqrt(sum(S_phi_f0(4:5))) / rms(phi))

delta_t = phi/(2*pi*1/nominal_interval);
figure;
plot(delta_t);
grid on;
xlabel('t / s');
ylabel('Time Fluctuation (TIE) / s')
title('Clock Difference in Time Domain')
