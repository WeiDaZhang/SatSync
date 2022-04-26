% TimeIntervalError from Phase Noise Power Law
function t_prime = PhaseNoisePowerLaw2TIE(N, nominal_interval, s_phi_f0)
% N is the length of the time series
% nominal_interval (s) is the ideal interval
% s_phi_f0 is a vector of s_phi (phase noise) for each f^beta slop intercepting
% f = 1Hz.
% beta = [-4,0] are random-walk-FM, flicker-FM, white-FM, flicker-PM, white-PM.
% size(s_phi_f0) = [5,1].

if(mod(N,2))
    N = N + 1;
end
mag = (1:N/2).^(-1).*randn(1,N/2);
pha = 2*pi*rand(1,N/2);
s_phi = [mag.*exp(1j*pha) conj(flip(mag.*exp(1j*pha)))];
figure;loglog(abs(s_phi))
% ifft(s_phi)
% 
% for (i=1;i<=N/2;i++) {
% 		mag = pow(i+1.0,-beta/2) * RandomGaussian(0.0,1.0); // Note to self a number of years later, why "i+1"
% 		pha = TWOPI * RandomUniform();
% 		real[i] = mag * cos(pha);
% 		imag[i] = mag * sin(pha);
% 		real[N-i] =  real[i];
% 		imag[N-i] = -imag[i];
% 	}
% 	imag[N/2] = 0;