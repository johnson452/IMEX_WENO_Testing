% Compute the moments:
function [n,u,T] = moments(f,app)

% Reconstruct f(v)
dir = "v";
poly_array_v_tilde = reconstruct(f,dir,app);

% Make the kernels (ans reversed polys)
kernel_n = [1,0,0];
kernel_v = [0,1,0];
kernel_v_sq = [0,0,1];


% Integrate to get the moments
n = quad_eval_mom(poly_array_v_tilde,kernel_n,app);
nu = quad_eval_mom(poly_array_v_tilde,kernel_v,app);
nu_sq = quad_eval_mom(poly_array_v_tilde,kernel_v_sq,app);

% Compute n, v, T (kb, m = 1)
m0 = app.m0;
kb = app.kb;
u = nu./n;
T = (m0./(kb*n)).*(nu_sq - n.*u.*u);

end