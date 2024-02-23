function [fp_im_half_limit,fm_ip_half_limit] = WENO_flux_limit(...
    fp_im_half,fm_ip_half,f_bar,fp_im_half_limit,fm_ip_half_limit,grid)
% 5th order WENO limiting using Section 2.4 of:
% REF: https://www.sciencedirect.com/science/article/pii/S0021999109007165

% Grab useful constants:
Nx = grid.Nx;

% Gauss-Labatto Quadrature Points (Must sum to 1 on this -1/2 to 1/2 interval)
w1 = 1/12;  % -1/2
w2 = 5/12;  % -0.447214/2
w3 = 5/12;  % +0.447214/2
w4 = 1/12;  % +1/2

% Compute the v_th location velocity grid
for v_index = 1:grid.Nv

    % Compute reconstruction Polynomial:
    for j = 1:Nx

        % (3.36a,b) of: https://epubs.siam.org/doi/epdf/10.1137/17M1144362
        xi_j = (f_bar(j,v_index) - fp_im_half(j,v_index)*w1 - fm_ip_half(j,v_index)*w4)/(w2 + w3);
        mj = min([fm_ip_half(j,v_index),fp_im_half(j,v_index),xi_j]);
        sigma_j = sigma(mj,f_bar(j,v_index));
        fp_im_half_limit(j,v_index) = sigma_j*( fp_im_half(j,v_index) - f_bar(j,v_index) ) + f_bar(j,v_index);
        fm_ip_half_limit(j,v_index) = sigma_j*( fm_ip_half(j,v_index) - f_bar(j,v_index) ) + f_bar(j,v_index);
    end
end
end



% Limit the poly nomial
function val = sigma(mj,f_bar_j)

% Compute sigma via 2.8
val = min([abs((f_bar_j)/(mj - f_bar_j)),1]);

end