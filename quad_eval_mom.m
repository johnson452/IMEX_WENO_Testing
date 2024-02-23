
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
        a = v(j)-dv/2;
        b = v(j)+dv/2;
        xj = v(j);

        % Grab the local polynomial
        poly_local = poly_array(:,i,j);

        % Iterate over quad-points
        for k = 1:N_quad

            % See: https://en.wikipedia.org/wiki/Gaussian_quadrature
            xi_norm = ((b-a)/2)*quad_x(k) + ((b+a)/2);
            kernel_val = sum(kernel.*[1,xi_norm,xi_norm^2]);
            mom(i) = mom(i) + ((b-a)/2)*w_quad(k)*poly_eval(poly_local,xj,xi_norm)*kernel_val;

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
