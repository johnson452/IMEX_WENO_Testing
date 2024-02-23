function eps = Knudsen_num(eps0,x)
% Knudsen Number

% Evaluate the knudsen number:
eps = eps0 + tanh(1 - 11*(x-1)) + tanh(1 + 11*(x-1));
end