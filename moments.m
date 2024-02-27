% Compute the moments:
function [n,u,T] = moments(f,app)



% Integrate to get the moments
if strcmp(app.grid.moments_type, "WENO_Reconstructed_fv") == 1

    % Reconstruct f(v)
    dir = "v";
    poly_array_v_tilde = reconstruct(f,dir,app);

    % Make the kernels (ans reversed polys)
    kernel_n = [1,0,0];
    kernel_v = [0,1,0];
    kernel_v_sq = [0,0,1];

    % Compute the polynomials
    n = quad_eval_mom(poly_array_v_tilde,kernel_n,app);
    nu = quad_eval_mom(poly_array_v_tilde,kernel_v,app);
    nu_sq = quad_eval_mom(poly_array_v_tilde,kernel_v_sq,app);
elseif strcmp(app.grid.moments_type, "Simple_No_Weno_reconst_fv") == 1
    [n,nu,nu_sq] = simple_moments(f,app);
end

% Compute n, v, T (kb, m = 1)
m0 = app.m0;
kb = app.kb;
u = nu./n;
T = (m0./(kb*n)).*(nu_sq - n.*u.*u);

end



% Simple moments where we do not polynomial reconstruct in velocity space
function [n,nu,nu_sq] = simple_moments(f,app)

%sz f
sz_f = size(f);

% Grab quantities
v = app.grid.v;
dv = app.grid.dv;

% Build moment arrays
n = zeros(sz_f(1),1);
nu = zeros(sz_f(1),1);
nu_sq = zeros(sz_f(1),1);

% Iterate over the whole grid
for i = 1:sz_f(1)
    n(i) = dv*sum(f(i,:));
    nu(i) = dv*sum(v.*f(i,:));
    nu_sq(i) = dv*sum(v.*v.*f(i,:));
end

end