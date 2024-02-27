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
%TEMP: ISOLATE G integration:
f_star = app.f; % - (dt/dx)*(app.F_jp_half - app.F_jm_half);

% (3.34b) Compute the moments (at ^{n + 1})
[n,u,T] = moments(f_star,app);
[n_func_tilde,u_func_tilde,T_func_tilde] = avg_to_func_moms(n,u,T,app);
f_Eq = equilibrium(n_func_tilde,u_func_tilde,T_func_tilde,app);

% Also save moments for diagnostics:
app.n_midpoint1 = n;
app.u_midpoint1 = u;
app.T_midpoint1 = T;


% Also save moments for diagnostics:
[n,u,T] = moments(squeeze(sum(f_Eq,1))/3,app);
app.n_Eq = n;
app.u_Eq = u;
app.T_Eq = T;

% (3.34c) Quadrature evaluation, need f_star_k(x), eps(c), M(x), f_Eq(x)
% via the reconstructed polynomials
f_star_func = reconstruct(f_star,"x",app);
%f_Eq_func = reconstruct(f_Eq,"x",app);
app.f = (1/dx)*quad_eval_int(f_star_func,f_Eq, app);

% Also save moments for diagnostics:
[n,u,T] = moments(app.f,app);
app.n = n;
app.u = u;
app.T = T;

%%% END IMEX SPP-Rk1 %%%


end


% Compute the quadrature evaluation
function [f_eval] = quad_eval_int(f_star_func,f_Eq,app)

% Grab constants from app
Nx = app.grid.Nx;
Nv = app.grid.Nv;
dx = app.grid.dx;
x = app.grid.x;
eps0 = app.eps0;

% 3-point Gauss-Legandre Quad: [-1 to 1]
N_quad = 3;
w_quad = [5/9,8/9,5/9];
xi_quad = [-sqrt(3/5),0,sqrt(3/5)];

% Build the equilibrium
f_eval = zeros(Nx,Nv);
dt = app.grid.dt;

% Iterate over the whole grid:
for i = 1:Nx

    % Bounds
    a = x(i)-dx/2;
    b = x(i)+dx/2;
    xj = x(i);

    for j = 1:Nv
        for k = 1:N_quad


            % See: https://en.wikipedia.org/wiki/Gaussian_quadrature
            xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);

            % Grab the values of f_star, f_Eq at quad points
            f_star = poly_eval(f_star_func(:,i,j),xj,xi_norm);
            %f_Eq = poly_eval(f_Eq_func(:,i,j),xj,xi_norm);

            % Grab the Knudsen number, compute integrand
            eps = Knudsen_num(eps0,xi_norm);
            integrand =  f_star/(1 + dt/eps) + ((dt/eps)/(1 + dt/eps))*f_Eq(k,i,j);

            % Equation 3.44 of
            % https://epubs.siam.org/doi/epdf/10.1137/17M1144362
            f_eval(i,j) = f_eval(i,j) + ((b-a)/2)*w_quad(k)*integrand;
        end
    end
end
end



% Compute the equilibrium everywhere
function f_Eq = equilibrium(n_func_tilde,u_func_tilde,T_func_tilde,app)

% Grab constants from app
Nx = app.grid.Nx;
Nv = app.grid.Nv;
v = app.grid.v;
dx = app.grid.dx;
x = app.grid.x;

% 3-point Gauss-Legandre Quad: [-1 to 1]
N_quad = 3;
w_quad = [5/9,8/9,5/9];
xi_quad = [-sqrt(3/5),0,sqrt(3/5)];

% Build the equilibrium
f_Eq = zeros(N_quad,Nx,Nv);

% Iterate over the whole grid:
for i = 1:Nx

    % Bounds
    a = x(i)-dx/2;
    b = x(i)+dx/2;
    xj = x(i);

    for k = 1:N_quad

        % Make the distribution at the quad points
        M_Eq = zeros(1,Nv);

        % See: https://en.wikipedia.org/wiki/Gaussian_quadrature
        xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);

        % Find n, u, T at the quad points
        n = poly_eval(n_func_tilde(:,i),xj,xi_norm);
        u = poly_eval(u_func_tilde(:,i),xj,xi_norm);
        T = poly_eval(T_func_tilde(:,i),xj,xi_norm);

        for j = 1:Nv

            % Equation 3.42 of
            % https://epubs.siam.org/doi/epdf/10.1137/17M1144362
            % Makes M^{n+1}/k
            M_Eq(j) = M_Eq(j) + maxwellian(n,u,T,v(j),app);
        end

        % Fix M_Eq
        M_Eq = fix_max(M_Eq,n,u,T,app);
        f_Eq(k,i,:) = M_Eq; %(1/2)*w_quad(k)*M_Eq;

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



% Compute fluxes - MAY NEED TO GENERALIZE FOR OTHER RK METHODS
function app = compute_fluxes(app)

% Direction:
dir = "x";

% Compute limited values: fp_im_half_limit, fm_ip_half_limit
[app.fp_im_half,app.fm_ip_half,app.f_bar] = WENO(app.f,app.fp_im_half,app.fm_ip_half,dir,app.grid);
[app.fp_im_half_limit,app.fm_ip_half_limit] = WENO_flux_limit(app.fp_im_half,app.fm_ip_half,app.f_bar,app.fp_im_half_limit,app.fm_ip_half_limit,app.grid);

% Compute fluxes via upwind 3.37
for i = 1:app.grid.Nx
    for j = 1:app.grid.Nv
        app.F_jp_half(i,j) = app.grid.v_plus(j)*app.fm_ip_half_limit(i,j) + app.grid.v_minus(j)*app.fp_im_half_limit(app.grid.R(i),j);
    end
end

% Periodic domain, build left
app.F_jm_half = app.F_jp_half(app.grid.L,:);

end
