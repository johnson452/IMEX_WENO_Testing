function eps = Knudsen_num(eps0,x)
% Knudsen Number

% Evaluate the knudsen number:
eps = eps0 + tanh(1 - 11*(x-1)) + tanh(1 + 11*(x-1));

% Diagnostic plot - matches fig 3 of J. Hu 18.
% x = linspace(0,2,100);
% eps = eps0 + tanh(1 - 11*(x-1)) + tanh(1 + 11*(x-1));
% plot(x,eps)
end