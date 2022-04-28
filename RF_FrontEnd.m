% RF front-end assumes signal has been sampled
% at f0 - ChipRate*nOverSample, and decimated downto ChipRate*nOverSample
% The front-end I&Q downconvert the input sigal at specified frequency

function IQBB = RF_FrontEnd(s_RF, ChipRate, nOverSample, f_IF_in)
N = length(s_RF);
t = 0 : 1/ChipRate/nOverSample : (N-1)/ChipRate/nOverSample;
f_IF = 1/4 * ChipRate * nOverSample;
if nargin > 3
    f_IF = f_IF_in;
end
IQBB = s_RF .* exp(1j*2*pi*f_IF*t);
f_IQBB = fft(IQBB);
IQBB = ifft([f_IQBB(1:floor(N/4)) zeros(1, N-2*floor(N/4)) f_IQBB(N-floor(N/4)+1:N)]);