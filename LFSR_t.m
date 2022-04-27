% Linear Feedback Shift Register by Time
function seq = LFSR_t(t, rate, state, tap)
% t (s) is array of time in seconds
% rate (Hz) is the symbol rate of LFSR
% state is the vector of initial state
% tap is the feedback polinomial, i.e. x^5 +x^2 + 1 == [5 2]

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
seq = zeros(1,N);
rounds = 0;
for idx = 1:N
    last_round = rounds;
    rounds = floor(t(idx)*rate);
    for idx_round = 1:(rounds - last_round)
        b = xor(state(tap(1)),state(tap(2)));
        if(long_tap)
            for idx_tap = 3:length(tap)
                b = xor(b, state(tap(idx_tap)));
            end
        end
        state = [b state(1:end - 1)];
    end
    seq(idx) = state(end);
end
    