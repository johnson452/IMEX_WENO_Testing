function app = SSP_IMEX_FIRST_ORDER(app)
% First order IMEX: (Eqns 3.33/3.34a-c):
% https://epubs.siam.org/doi/epdf/10.1137/17M1144362

% Gran useful quantities
grid = app.grid;
dx = grid.dx;
dt = grid.dt;

%%% START IMEX SPP-Rk1 %%%
% (3.34a) Compute f_{j,k}^star:
app = compute_fluxes(app);
f_star = app.f - (dt/dx)*(app.F_jp_half - app.F_jm_half);

% (3.34b) Compute the moments (at ^{n + 1})
[n,u,T] = moments(f_star,app);
[n_func_tilde,u_func_tilde,T_func_tilde] = avg_to_func_moms(n,u,T,app);
f_Eq = equilibrium(n_func_tilde,u_func_tilde,T_func_tilde,app);

% (3.34c) Quadrature evaluation, need f_star_k(x), eps(c), M(x), f_Eq(x)
% via the reconstructed polynomials
app.f = (1/dx)*quad_eval(   f_star_func/(1 + dt/eps_func) + ...
    ((dt/eps_func)/(1 + dt/eps_func))*f_Eq_func   );
%%% END IMEX SPP-Rk1 %%%


end



% Compute the equilibiurm everywhere
function f_Eq = equilibrium(n_func_tilde,u_func_tilde,T_func_tilde,app)

% Grab constants from app
Nx = app.grid.Nx;
Nv = app.grid.Nv;
v = app.grid.v;
dx = app.grid.dx;
x = app.grid.x;

% 3-point Gauss-Legandre Quad: [-1/2 to 1/2]
N_quad = 3;
w_quad = [5/9,8/9,5/9]/2;
xi_quad = [-sqrt(3/5),0,sqrt(3/5)]/2;

% Build the equilibrium
f_Eq = zeros(Nx,Nv);

% Iterate over the whole grid:
for i = 1:Nx
    for j = 1:Nv
        for k = 1:N_quad

            % Find n, u, T at the quad points
            n = poly_eval(n_func_tilde(:,i),x(i),xi_quad(k)*dx);
            u = poly_eval(u_func_tilde(:,i),x(i),xi_quad(k)*dx);
            T = poly_eval(T_func_tilde(:,i),x(i),xi_quad(k)*dx);

            % Equation 3.42 of
            % https://epubs.siam.org/doi/epdf/10.1137/17M1144362 
            % Makes M^{n+1}/k
            f_Eq(i,j) = f_Eq(i,j) + w_quad(k)*maxwellian(n,u,T,v(j),app);
        end
    end
end


end


% Function, averges to functions
function [n_func_tilde,u_func_tilde,T_func_tilde] = avg_to_func_moms(n,u,T,app)

% Construct the functions:
dir = "x";
n_func_tilde = reconstruct(n,dir,app);
u_func_tilde = reconstruct(u,dir,app);
T_func_tilde = reconstruct(T,dir,app);


end


% Reconstruct
function [poly_array_tilde] = reconstruct(f,dir,app)

% WENO reconstruct -> gives f_{bar}, fp_{kp_half}, fm_{km_half}
fp_km_half = zeros(size(f));
fm_kp_half = zeros(size(f));
[fp_km_half,fm_kp_half,f_bar] = WENO(f,fp_km_half,fm_kp_half,dir,app.grid);
poly_array_tilde = WENO_poly_build_and_limit(fp_km_half,fm_kp_half,f_bar,dir,app.grid);

end


% Compute the moments:
function [n,u,T] = moments(f,app)

% Reconstruct f(v)
dir = "v";
poly_array_v_tilde = reconstruct(f,dir,app);

% Make the kernels (ans reversed polys)
kernel_n = 1;
kernel_v = [0,1];
kernel_v_sq = [0,0,1];


% Integrate to get the moments
n = quad_eval_mom(poly_array_v_tilde,kernel_n,app);
nv = quad_eval_mom(poly_array_v_tilde,kernel_v,app);
nv_sq = quad_eval_mom(poly_array_v_tilde,kernel_v_sq,app);

% Compute n, v, T (kb, m = 1)
m0 = app.m0;
kb = app.kb;
u = nv./n;
T = (m0./(kb*n)).*(nv_sq - n.*u.*u);

end


% Compute the quadrature evaluation
function [mom] = quad_eval_mom(poly_array,kernel,app)

% Fix:
% Grab useful constants:
Nx = app.grid.Nx;
dv = app.grid.dv;
Nv = app.grid.Nv;
v = app.grid.v;


% Setup the quad:
% So Gauss-Legandre is exact to poly order 2n-1 (since kernels can take x^4
% -> x^6
N_quad = 4;
% Domain -1 to 1
quad_x = [-sqrt(3/7 - (2/7)*sqrt(6/5)), -sqrt(3/7 + (2/7)*sqrt(6/5)), sqrt(3/7 + (2/7)*sqrt(6/5)), sqrt(3/7 - (2/7)*sqrt(6/5))];
w_quad = [(18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36,(18+sqrt(30))/36];

% Biuld the moment array
mom = zeros(Nx,1);

% Iterate over the x-grid
for i = 1:Nx

    % Iterate over the velocity grid:
    for j = 1:Nv

        % Bounds
        a = v(i)-dv/2;
        b = v(i)-dv/2;
        xj = v(i);

        % Grab the local polynomial
        poly_local = poly_array(:,i,j);

        % Multiply the kernel by the polynomial
        poly_local = poly_multiply(poly_local,kernel);

        % Iterate over quad-points
        for k = 1:N_quad

            % See: https://en.wikipedia.org/wiki/Gaussian_quadrature
            xi_norm = ((b-a)/2)*quad_x(k) + ((b+a)/2);
            mom(i) = mom(i) + ((b-a)/2)*w_quad(k)*poly_eval(poly_local,xj,xi_norm);

        end
    end
end
end


% Multiply the kernel:
function [poly_res] = poly_multiply(p1,p2)

% Multiply the kernel
sz_poly1 = max(size(p1));
sz_poly2 = max(size(p2));
poly_res = zeros(sz_poly2+sz_poly1 - 1,1);

% Iterate over all combinations
for i = 1:sz_poly1
    for j = 1:sz_poly2
        poly_res(i+j-1) = poly_res(i+j-1) + p1(i)*p2(j);
    end
end

end


% Compute fluxes - MAY NEED TO GENERALIZE FOR OTHER RK METHODS
function app = compute_fluxes(app)

% Compute limited values: fp_im_half_limit, fm_ip_half_limit
[app.fp_im_half,app.fm_ip_half,app.f_bar] = WENO(app.f,app.fp_im_half,app.fm_ip_half,app.grid);
[app.fp_im_half_limit,app.fm_ip_half_limit] = WENO_flux_limit(app.fp_im_half,app.fm_ip_half,app.f_bar,app.fp_im_half_limit,app.fm_ip_half_limit,app.grid);

% Compute fluxes via upwind 3.37
app.F_jp_half = app.grid.v_plus.*app.fm_ip_half_limit + app.grid.v_minus.*app.fp_im_half_limit(app.grid.R,:);

% Periodic domain, build left
app.F_jm_half = app.F_jp_half(app.grid.L,:);

end