% Rel_diff
function res = rel_diff(a,b)

% Compute the relative difference between a and b
res = abs(a - b)/abs((a+b)/2);

end