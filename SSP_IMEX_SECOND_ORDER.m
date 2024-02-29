function app = SSP_IMEX_SECOND_ORDER(app)
% Second order IMEX: From Dingyun's writeup Eq. 9

% Gran useful quantities
grid = app.grid;
dx = grid.dx;
dt = grid.dt;

% Butcher coefficients (2.6.1 J. Hu 2018)
a_tilde_21 = 0.73695027152854;
a_tilde_31 = 0.32152816910844;
a_tilde_32 = 0.67847183089156;
a_11 = 0.62863517121833;
a_21 = 0.24310046553707;
a_31 = 0.48036510509894;
a_22 = 0.19593925696632;
a_32 = 0.074643281386981;
a_33 = 0.44499161351408;
alpha = 0.27973737915215;


%%% START IMEX SPP-Rk2,  %%%
 %%% STAGE 1:
% (9a - Sub-Stage a) Compute f_{j,k}^star:
f_star = app.f; % - (dt/dx)*(app.F_jp_half - app.F_jm_half);

% (9a - Sub-Stage b) Compute the moments (at ^{n + 1})
[n1,u1,T1] = moments(f_star,app);
[n1_func_tilde,u1_func_tilde,T1_func_tilde] = avg_to_func_moms(n1,u1,T1,app);
f1_Eq = equilibrium(n1_func_tilde,u1_func_tilde,T1_func_tilde,app);

% Also save moments for diagnostics:
app.n_midpoint1 = n1; app.u_midpoint1 = u1; app.T_midpoint1 = T1;

% (9a - Sub-Stage c) Quadrature evaluation, need f_star_k(x), eps(c), M(x), f_Eq(x)
% via the reconstructed polynomials
f_star_func = reconstruct(f_star,"x",app);
f1 = (1/dx)*quad_stg_k(f_star_func,f1_Eq,a_11,app);


%%% STAGE 2:
% (9b - Sub-Stage a) Compute f_{j,k}^star:
f1_func = reconstruct(f1,"x",app);
[F1_jp_half,F1_jm_half] = compute_fluxes(f1,app);
f_star_star = app.f - (dt/dx)*a_tilde_21*(F1_jp_half - F1_jm_half) + ...
    (dt/dx)*a_21*quad_Eq(f1_func,f1_Eq,app);

% (9b - Sub-Stage b) Compute the moments (at ^{n + 1})
[n2,u2,T2] = moments(f_star_star,app);
[n2_func_tilde,u2_func_tilde,T2_func_tilde] = avg_to_func_moms(n2,u2,T2,app);
f_Eq2 = equilibrium(n2_func_tilde,u2_func_tilde,T2_func_tilde,app);

% Also save moments for diagnostics:
app.n_midpoint2 = n2; app.u_midpoint2 = u2; app.T_midpoint2 = T2;

% (9b - Sub-Stage c) Quadrature evaluation, need f_star_k(x), eps(c), M(x), f_Eq(x)
% via the reconstructed polynomials
f_star_star_func = reconstruct(f_star_star,"x",app);
f2 = (1/dx)*quad_stg_k(f_star_star_func,f_Eq2,a_22,app);


%%% STAGE 3:
% (9c - Sub-Stage a) Compute f_{j,k}^star:
f2_func = reconstruct(f2,"x",app);
[F2_jp_half,F2_jm_half] = compute_fluxes(f2,app);
f_star_star_star = app.f - (dt/dx)*a_tilde_31*(F1_jp_half - F1_jm_half) + ...
    (dt/dx)*a_31*quad_Eq(f1_func,f1_Eq,app) ...
    - (dt/dx)*a_tilde_32*(F2_jp_half - F2_jm_half) + ...
    (dt/dx)*a_32*quad_Eq(f2_func,f_Eq2,app);

% (9b - Sub-Stage b) Compute the moments (at ^{n + 1})
[n3,u3,T3] = moments(f_star_star_star,app);
[n3_func_tilde,u3_func_tilde,T3_func_tilde] = avg_to_func_moms(n3,u3,T3,app);
f_Eq3 = equilibrium(n3_func_tilde,u3_func_tilde,T3_func_tilde,app);

% Also save moments for diagnostics:
app.n_midpoint3 = n3; app.u_midpoint3 = u3; app.T_midpoint3 = T3;

% (9b - Sub-Stage c) Quadrature evaluation, need f_star_k(x), eps(c), M(x), f_Eq(x)
% via the reconstructed polynomials
f_star_star_star_func = reconstruct(f_star_star_star,"x",app);
f3 = (1/dx)*quad_stg_k(f_star_star_star_func,f_Eq3,a_33,app);


% Finally compute the last stage:
f3_func = reconstruct(f3,"x",app);
app.f = (1/dx)*quad_stg_final(f3_func,f_Eq3,alpha,app);
%%% END IMEX SPP-Rk2 %%%


% Also save moments for diagnostics:
[n,u,T] = moments(app.f,app);
app.n = n;
app.u = u;
app.T = T;

end


% Compute the quadrature evaluation
function [f_eval] = quad_Eq(fk_func,f_Eq,app)

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
%dt = app.grid.dt;

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
            fk = poly_eval(fk_func(:,i,j),xj,xi_norm);
            %f_Eq = poly_eval(f_Eq_func(:,i,j),xj,xi_norm);

            % Grab the Knudsen number, compute integrand
            eps = Knudsen_num(eps0,xi_norm);
            integrand =  (1/eps) * ( f_Eq(k,i,j) - fk );

            % Equation 3.44 of
            % https://epubs.siam.org/doi/epdf/10.1137/17M1144362
            f_eval(i,j) = f_eval(i,j) + ((b-a)/2)*w_quad(k)*integrand;
        end
    end
end
end


% Compute the quadrature evaluation
function [f_eval] = quad_stg_final(f_nu,f_Eq,alpha,app)

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
            f_nu_val = poly_eval(f_nu(:,i,j),xj,xi_norm);
            %f_Eq = poly_eval(f_Eq_func(:,i,j),xj,xi_norm);

            % Grab the Knudsen number, compute integrand
            eps = Knudsen_num(eps0,xi_norm);
            integrand =  (1/(1 + alpha*dt*dt/(eps^2) )) * ( f_nu_val + alpha*((dt/eps)^2)*f_Eq(k,i,j ));

            % Equation 3.44 of
            % https://epubs.siam.org/doi/epdf/10.1137/17M1144362
            f_eval(i,j) = f_eval(i,j) + ((b-a)/2)*w_quad(k)*integrand;
        end
    end
end
end



% Compute the quadrature evaluation
function [f_eval] = quad_stg_k(f_star_func,f_Eq,aii,app)

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
            integrand =  (1/(1 + aii*dt/eps)) * ( f_star + aii*((dt/eps))*f_Eq(k,i,j ));

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
%w_quad = [5/9,8/9,5/9];
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
function [F_jp_half,F_jm_half] = compute_fluxes(f,app)

% Direction:
dir = "x";
F_jp_half = zeros(size(f));
fp_im_half = zeros(size(f));
fm_ip_half = zeros(size(f));
fp_im_half_limit  = zeros(size(f)); 
fm_ip_half_limit = zeros(size(f));

% Compute limited values: fp_im_half_limit, fm_ip_half_limit
[fp_im_half,fm_ip_half,f_bar] = WENO(f,fp_im_half,fm_ip_half,dir,app.grid);
[fp_im_half_limit,fm_ip_half_limit] = WENO_flux_limit(fp_im_half,fm_ip_half,f_bar,fp_im_half_limit,fm_ip_half_limit,app.grid);

% Compute fluxes via upwind 3.37
for i = 1:app.grid.Nx
    for j = 1:app.grid.Nv
        F_jp_half(i,j) = app.grid.v_plus(j)*fm_ip_half_limit(i,j) + app.grid.v_minus(j)*fp_im_half_limit(app.grid.R(i),j);
    end
end

% Periodic domain, build left
F_jm_half = F_jp_half(app.grid.L,:);

end
