function [poly_array_tilde] = WENO_poly_build_and_limit(...
    fp_im_half,fm_ip_half,f_bar,dir,grid)
% 5th order WENO limiting using Section 2.4 of:
% REF: https://www.sciencedirect.com/science/article/pii/S0021999109007165

% Grab size of the array
sz_f = size(f_bar);

% Grab useful constants:
Nx = sz_f(1);
Nv = sz_f(2);
dx = grid.dx;

% If the direction is v:
if strcmp(dir,"v") == 1
    Nx = sz_f(2);
    dx = grid.dv;
    Nv = sz_f(1);
    fp_im_half = transpose(fp_im_half);
    fm_ip_half = transpose(fm_ip_half);
    f_bar = transpose(f_bar);
end

% Gauss-Labatto Quadrature Points (Must sum to 1 on this -1/2 to 1/2 interval)
w1 = 1/12;  % -1/2
w2 = 5/12;  % -0.447214/2
w3 = 5/12;  % +0.447214/2
w4 = 1/12;  % +1/2

% Decide the quad. points
%Sj = [-1/2,-1/sqrt(5),1/sqrt(5),1/2]*dx;

% Build the poly array: [a0, ... , a5]
poly_array_tilde = zeros(5,Nx,Nv);

% Compute the v_th location velocity grid
for v_index = 1:Nv

    % Compute reconstruction Polynomial:
    for j = 1:Nx

        % Periodic domains:
        jm = mapindex(j-1,Nx);
        jp = mapindex(j+1,Nx);

        % Compute poly coeff
        a0 = (1/192)*(f_bar(jm) + 298*f_bar(j) + f_bar(jp) - 54*(fp_im_half(j) + fm_ip_half(j)));
        a1 = (1/(8*dx))*(f_bar(jm) - f_bar(jp) - 10*(fp_im_half(j) - fm_ip_half(j)));
        a2 = (1/(8*dx^2))*( -(f_bar(jm) + 58*f_bar(j) + f_bar(jp)) + 30*(fp_im_half(j) + fm_ip_half(j)));
        a3 = (1/(8*dx^3))*(f_bar(jp) - f_bar(jm) +2*(fp_im_half(j) - fm_ip_half(j)));
        a4 = (1/(12*dx^4))*( (5*f_bar(jm) + 50*f_bar(j) + 5*f_bar(jp)) - 30*(fp_im_half(j) + fm_ip_half(j)));

        % Poly eval:
        poly = [a0,a1,a2,a3,a4];

        % (3.35a,b) of: https://epubs.siam.org/doi/epdf/10.1137/17M1144362
        xi_j = (f_bar(j,v_index) - fp_im_half(j,v_index)*w1 - fm_ip_half(j,v_index)*w4)/(w2 + w3);
        mj = min([fm_ip_half(j,v_index),fp_im_half(j,v_index),xi_j]); % equiv to eval_poly, see PP before
        sigma_j = sigma(mj,f_bar(j,v_index));

        %Limit the polynomial
        poly_tilde = sigma_j*(poly) - sigma_j*[1,0,0,0,0]*f_bar(j,v_index) + [1,0,0,0,0]*f_bar(j,v_index);

        % Save the polynomial
        poly_array_tilde(:,j,v_index) = poly_tilde;
    end
end

% Permute the last two arrays:
% If the direction is v:
if strcmp(dir,"v") == 1
    poly_array_tilde = permute(poly_array_tilde, [1, 3, 2]);
end

end


% Limit the poly nomial
function val = sigma(mj,f_bar_j)

% Compute sigma via 2.8
val = min([abs((f_bar_j)/(mj - f_bar_j)),1]);

end


% Map the the grid indices
function [x] = mapindex(x,Nx)
for i = 1:max(size(x))
    if (x(i) > Nx)
        x(i) = x(i) - Nx;
    end
    if (x(i) < 1)
        x(i) = x(i) + Nx;
    end
end
end