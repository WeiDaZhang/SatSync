% Linear Feedback Shift Register by Time
function seq = LFSR_t(t, rate, state, tap, chip_length)
% t (s) is array of time in seconds
% rate (Hz) is the symbol rate of LFSR
% state is the vector of initial state
% tap is the feedback polinomial, i.e. x^5 +x^2 + 1 == [5 2]
% chip_length is the length of the chip

if(max(tap) > length(state))
    MSB = max(tap);
    if(size(state,1) == 1)
        state = [zeros(1,MSB - length(state)) state];
    else
        state = [zeros(MSB - length(state),1); state]';
    end
else
    MSB = length(state);
end

if(tap(end) == 0)
    tap = tap(1:end-1);
end

if(isempty(find(state)))
    disp('Initial state changed to [zeros 1]')
end

N = length(t);
long_tap = length(tap) > 2;     % extra xor needed if tap number > 2
b = false;
chip = zeros(1, chip_length);
seq = zeros(1, N);
for idx_chip = 1:chip_length
    b = xor(state(tap(1)),state(tap(2)));
    if(long_tap)
        for idx_tap = 3:length(tap)
            b = xor(b, state(tap(idx_tap)));
        end
    end
    state = [b state(1:end - 1)];
    chip(idx_chip) = state(end);
end

for idx = 1:N
    seq(idx) = chip(mod(floor(t(idx)*rate), chip_length) + 1)*2-1;
end
