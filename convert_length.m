function [vector_of_extend] = convert_length(vector_of_length)
n = length(vector_of_length);
idx = find(vector_of_length ~= 1);
oneidx = find(vector_of_length == 1);
vector_of_extend = ones(1, n);
if(isempty(idx))
    return
elseif(max(vector_of_length(idx)) == min(vector_of_length(idx)))
    vector_of_extend(oneidx) = max(vector_of_length(idx));
    return
else
    vector_of_extend = zeros(1,n);
end
