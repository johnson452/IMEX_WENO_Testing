function [fp_im_half,fm_ip_half,f_bar_j] = WENO(f_full,fp_im_half,fm_ip_half,dir,grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes f (cell averages) and the grid, and returns
% Cell: |i-1/2     i     i+1/2|
%       |fp       fi        fm|

% f_full is expected to be a phase space/phase-space flux: f(x,v) or
% v*f(x,v)

% Only works for 5th order!

%From: https://www.sciencedirect.com/science/article/pii/S0021999109007165
% and: https://www3.nd.edu/~zxu2/acms60790S13/Shu-WENO-notes.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab certain quantities
k = grid.WENO_order;
x = grid.x;
dx = grid.dx;
sz_f = size(f_full);
Nx = sz_f(1);
Nv = sz_f(2);
f_bar_j = zeros(sz_f);

% 3-point Gauss-Legandre Quad: [-1 to 1]
N_quad = 3;
w_quad = [5/9,8/9,5/9];
xi_quad = [-sqrt(3/5),0,sqrt(3/5)];

% Make sure we aren't calling this by accident
if strcmp(grid.moments_type,"Simple_No_Weno_reconst_fv") == 1 && strcmp(dir,"v") ==1
    printf("Incorrect call to WENO to reconstruct v!\n")
end


% Transpose the direction if dir = "v"
if dir == "v"
    x = grid.v;
    dx = grid.dv;
    Nx = sz_f(2);
    Nv = sz_f(1);
    f_full = transpose(f_full);
    f_bar_j = transpose(f_bar_j);
    fp_im_half = transpose(fp_im_half);
    fm_ip_half = transpose(fm_ip_half);
end

% Compute the v_th location velocity grid
for v_index = 1:Nv

    % Compute a local f
    f = f_full(:,v_index);

    % Compute the ith location spatially
    for i = 1:Nx

        %1. Obtain the k reconstructed values v(r) , of k-th order accuracy, in (2.51),
        %based on the stencils (2.50), for r = 0,...,k − 1;
        %Also obtain the k reconstructed values v(r) , of k-th order accuracy,
        %using (2.10), again based on i−1 the stencils (2.50), for r = 0, ..., k − 1;

        % Create a vector for f^r_p, f^r_m
        f_r_ip_half = zeros(1,k);
        f_r_im_half = zeros(1,k);
        Beta_r = zeros(1,k);
        alpha_r = zeros(1,k);
        alpha_r_tilde = zeros(1,k);
        wr = zeros(1,k);
        wr_tilde = zeros(1,k);

        for r = 0:k - 1

            % Construct the polynomial
            pre_xr_index = linspace(i-r,i-r+k-1,k);
            xr_index = mapindex(dir,pre_xr_index,Nx);
            yr_poly = f(xr_index);
            ujm = yr_poly(1);
            uj = yr_poly(2);
            ujp = yr_poly(3);

            %pr = polyfit(xr_poly,yr_poly,n); %2.50, constrained by 2.9
            a2 = (ujm - 2*uj + ujp)/(2*dx^2);
            a1 = -(ujm - ujp)/(2*dx);
            a0 = (13*uj)/12 - ujm/24 - ujp/24;
            pr = [a0,a1,a2];


            % Checkpoly
            %check_poly_second_order(pr,uj,ujm,ujp,x(i),dx)
            

            %             % Compute v_bar_j
            %             for j = i-r:i+(k - 1 - r) % s = k - 1 - r (2.9)
            %
            %                 % Set the initial value to zero
            %                 j_index = mapindex(dir,j,Nx);
            %                 f_bar_j(j_index,v_index) = 0;
            %
            %                 %Iterate over quad points
            %                 for k = 1:N_quad
            %
            %                     % Compute the mean value
            %                     xj = x_ghost(Nx,j_index,j,dir,dx,x);
            %                     a = xj - dx/2;
            %                     b = xj + dx/2;
            %                     xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);
            %                     f_bar_j(j_index,v_index) = f_bar_j(j_index,v_index) + ...
            %                         (1/dx)*((b-a)/2)*w_quad(k)*poly_eval(pr,xj,xi_norm);
            %
            %
            %                 end
            %
            %
            %                 % Check diff
            %                 if rel_diff(f_bar_j(j_index,v_index),f(j_index)) > 1e-8
            %                     fprintf("WENO fails to build poly!\n");
            %                 end
            %             end

            % Since we've shown equiv
            f_bar_j(xr_index,v_index) = f(xr_index);



            % (2.11) v_bar
            f_ip_half = 0;
            f_im_half = 0;
            for j = 0:k-1
                f_ip_half = f_ip_half + crj(k,r,j)*f_bar_j(mapindex(dir,i - r + j,Nx),v_index);
                f_im_half = f_im_half + crj_tilde(k,r,j)*f_bar_j(mapindex(dir,i - r + j,Nx),v_index);
            end

            % Save the output (r + 1) for matlab
            f_r_ip_half(r + 1) = f_ip_half;
            f_r_im_half(r + 1) = f_im_half;


            %3. Find the smooth indicators βr in (2.61), for all r = 0, k − 1.
            %Explicit formulae for k = 2 and k = 3 are given in (2.62) and (2.63) respectively.

            % Compute Beta_r, takes care of k = 1 case, sets Beta+_r(1) = 0
            % by default
            %             for l = 1:k-1
            %
            %                 % Compute 2.61
            %                 dp = polyder(pr,l);
            %                 p_squared = conv(dp, dp);
            %                 q = polyint(p_squared);
            %
            %                 % Save Beta_r to an array
            %                 a = x(i)-dx/2;
            %                 b = x(i)+dx/2;
            %                 Beta_r(r + 1) = Beta_r(r + 1) + (dx^(2*l-1))*diff(polyval(q,[a b]));
            %
            %             end

            %Compute beta_analytic
            Beta_r(r + 1) = beta_analytic_calc(i,f_bar_j(:,v_index),k,r,dir,Nx);

            % Check Beta: Checks out
            %check_beta(i,f_bar_j(:,v_index),k,Beta_r,r,Nx);

            %2. Find the constants d and d ̃ , such that (2.54) and (see nd.edu, pg 20)

            % Compute alpha_r
            eps = 1e-6;
            alpha_r(r+1) = dr(r,k)/((Beta_r(r + 1) + eps)^2);
            alpha_r_tilde(r+1) = dr_tilde(r,k)/((Beta_r(r + 1) + eps)^2);

        end %end r

        %4. Form the weights ωr and ω r using (2.58)-(2.59) and (see nd.edu, pg 20)

        for r = 0:k - 1

            % Compute the weights
            denom = 0;
            for l = 1:k
                denom = denom + alpha_r(l);
            end
            wr(r+1) = alpha_r(r+1)/denom;

            % Compute the weights - tilde
            denom_tilde = 0;
            for l = 1:k
                denom_tilde = denom_tilde + alpha_r_tilde(l);
            end
            wr_tilde(r+1) = alpha_r_tilde(r+1)/denom_tilde;

        end %end r


        %5. Find the (2k − 1)-th order reconstruction (2.64)
        fp_im_half(i,v_index) = 0;
        fm_ip_half(i,v_index) = 0;
        %check_wr(wr);
        %check_wr(wr_tilde);
        for r = 0:k-1
            fm_ip_half(i,v_index) = fm_ip_half(i,v_index) + wr(r+1)*f_r_ip_half(r + 1);
            fp_im_half(i,v_index) = fp_im_half(i,v_index) + wr_tilde(r+1)*f_r_im_half(r + 1);
        end

    end % End iteration over x

end % End iteration over v

% Transpose the direction if dir = "v"
if dir == "v"
    f_bar_j = transpose(f_bar_j);
    fp_im_half = transpose(fp_im_half);
    fm_ip_half = transpose(fm_ip_half);
end

end % End Function


% Grab the ghost cells:
function [x_vals] = x_ghost(Nx,fixed_indices,indices,dir,dx,x)

% Build ghost cells
x_vals = max(size(indices));
if strcmp(dir,"x") == 1
    x_vals = x(fixed_indices);
elseif (strcmp(dir,"v") == 1)
    for i = 1:max(size(indices))
        indx = indices(i);
        if indx > Nx
            x_vals(i) = x(Nx) + (indx - Nx)*dx;
        elseif indx < 1
            x_vals(i) = x(1) + (indx - 1 )*dx;
        else
            x_vals(i) = x(indx);
        end
    end
end
end


% % Function: Unit Test on Beta
% function check_beta(i,f_bar_j,k,Beta_r,r,dir,Nx)
%
% analytic = beta_analytic_calc(i,f_bar_j,k,r,dir,Nx);
%
% if (Beta_r(r + 1)-analytic)/((Beta_r(r + 1)+analytic)/2) > 1e-8 && ...
%         abs((Beta_r(r + 1)+analytic)/2) > 1e-8
%     fprintf("Beta fails to reproduce analytic results (k = %d)\n",k);
% end
%
%
% % Extra diagnostic:
% %fprintf("(2.61): %1.16e, (analytic): %1.16e, (Diff): %1.16e\n", ...
% %Beta_r(r + 1), analytic,Beta_r(r + 1)-analytic);
% end

% % Function: Unit test on Wr
% function check_wr(x)
%
% % Check the property
% if abs(sum(x) -1) > 1e-8
%     fprintf("wr ill formed!\n")
% end
%
% end

function analytic = beta_analytic_calc(i,f_bar_j,k,r,dir,Nx)
if k == 1
    analytic = 1;
elseif k == 2
    if r == 0
        analytic = (f_bar_j(mapindex(dir,i+1,Nx))- f_bar_j(i))^2;
    elseif r == 1
        analytic = (f_bar_j(i) - f_bar_j(mapindex(dir,i-1,Nx)))^2;
    end
elseif k == 3
    if r == 0
        analytic = (13/12)*(f_bar_j(i) - 2*f_bar_j(mapindex(dir,i+1,Nx)) + f_bar_j(mapindex(dir,i+2,Nx)))^2 +...
            (1/4)*(3*f_bar_j(i) - 4*f_bar_j(mapindex(dir,i+1,Nx)) + f_bar_j(mapindex(dir,i+2,Nx)))^2;
    elseif r == 1
        analytic = (13/12)*(f_bar_j(mapindex(dir,i-1,Nx)) - 2*f_bar_j(i) + f_bar_j(mapindex(dir,i+1,Nx)))^2 +...
            (1/4)*(f_bar_j(mapindex(dir,i-1,Nx)) - f_bar_j(mapindex(dir,i+1,Nx)))^2;
    elseif r == 2
        analytic = (13/12)*(f_bar_j(mapindex(dir,i-2,Nx)) - 2*f_bar_j(mapindex(dir,i-1,Nx)) + f_bar_j(i))^2 +...
            (1/4)*(f_bar_j(mapindex(dir,i-2,Nx)) - 4*f_bar_j(mapindex(dir,i-1,Nx)) + 3*f_bar_j(i))^2;
    end
end
end



% dr_tilde table
function [dr_tilde_val] = dr_tilde(r,k)
dr_tilde_val = dr(k-1-r,k);

end


% dr table
function [dr_val] = dr(r,k)

%Cases
if k == 1
    if r == 0
        dr_val = 1;
    end
elseif k == 2
    if r == 0
        dr_val = 2/3;
    end
    if r == 1
        dr_val = 1/3;
    end
elseif k == 3
    if r == 0
        dr_val = 3/10;
    end
    if r == 1
        dr_val = 3/5;
    end
    if r == 2
        dr_val = 1/10;
    end
else
    fprintf("dr is not defined for these values\n");
end

end


% crj_tilde table
function val = crj_tilde(k,r,j)
val = crj(k,r-1,j);
end


% crj table
function val = crj(k,r,j)

% From: nd.edu (table in section 2)
rm = r + 2;
jm = j + 1;

% k-case:
if k == 1

    % crj_matrix
    crj_matrix = [1;...
        1];

    % Cases of r and j
    if (r >= -1 && r <= k-1) && (j >= 0 && j<= k-1)
        val = crj_matrix(rm,jm);
    else
        fprintf("crj value doesn't exist! (%d) \n",k);
    end


elseif k == 2

    % crj_matrix
    crj_matrix = [3/2,-1/2;...
        1/2,1/2;...
        -1/2,3/2];

    % Cases of r and j
    if (r >= -1 && r <= k-1) && (j >= 0 && j<= k-1)
        val = crj_matrix(rm,jm);
    else
        fprintf("crj value doesn't exist! (%d) \n",k);
    end

elseif k == 3

    % crj_matrix
    crj_matrix = [11/6,-7/6, 1/3;...
        1/3,5/6,-1/6;...
        -1/6,5/6,1/3;...
        1/3,-7/6,11/6];

    % Cases of r and j
    if (r >= -1 && r <= k-1) && (j >= 0 && j<= k-1)
        val = crj_matrix(rm,jm);
    else
        fprintf("crj value doesn't exist! (%d) \n",k);
    end

else
    fprintf("No crj specified for this {k,r,j}\n");
end

%Verify matrix with: Done!
%should_be_zero = max(crj_matrix - flip(flip(crj_matrix,2),1));
%check_crj(should_be_zero);


end


% function check_crj(x)
%
% if max(x) ~= 0 || min(x) ~= 0
%     fprintf("Issue in crj_matrices\n");
% end
%
% end



% Check the polynomial fourth order matches the condtions
function check_poly_second_order(poly,f_bar,fL_bar,fR_bar,xj,dx)

% Does the edge evaluations match?
x_jm_half = xj - dx/2;
x_jp_half = xj + dx/2;

% 3-point Gauss-Legandre Quad: [-1 to 1]
N_quad = 3;
w_quad = [5/9,8/9,5/9];
xi_quad = [-sqrt(3/5),0,sqrt(3/5)];

% Check the means:
a = x_jm_half;
b = x_jp_half;
int_I = 0;
for k = 1:N_quad
    xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);
    int_I = int_I + (1/dx)*((b-a)/2)*w_quad(k)*poly_eval(poly,xj,xi_norm);
end
err_int = rel_diff(int_I,f_bar);

% Check the means:
a = x_jm_half-dx;
b = x_jp_half-dx;
int_L = 0;
for k = 1:N_quad
    xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);
    int_L = int_L + (1/dx)*((b-a)/2)*w_quad(k)*poly_eval(poly,xj,xi_norm);
end
err_L_int = rel_diff(int_L,fL_bar);

% Check the means:
a = x_jm_half+dx;
b = x_jp_half+dx;
int_R = 0;
for k = 1:N_quad
    xi_norm = ((b-a)/2)*xi_quad(k) + ((b+a)/2);
    int_R = int_R + (1/dx)*((b-a)/2)*w_quad(k)*poly_eval(poly,xj,xi_norm);
end
err_R_int = rel_diff(int_R,fR_bar);

tol = 1e-8;

if err_int > tol || err_L_int > tol || err_R_int > tol

    fprintf("(POLY 2nd-O MEAN) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_I,f_bar,err_int);
    fprintf("(POLY 2nd-O MEAN L) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_L,fL_bar,err_L_int);
    fprintf("(POLY 2nd-O MEAN R) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_R,fR_bar,err_R_int);
    fprintf("\n");

    clf()
    x_lower = xj - 3*dx/2;
    x_upper = xj + 3*dx/2;
    x_span = linspace(x_lower,x_upper);
    y_span = poly_eval(poly,xj,x_span);
    plot(x_span,y_span,"red")
    hold on
    plot([xj - 3*dx/2,xj - dx/2],[fL_bar,fL_bar],"black")
    hold on
    plot([xj - dx/2,xj + dx/2],[f_bar,f_bar],"black")
    hold on
    plot([xj + dx/2,xj + 3*dx/2],[fR_bar,fR_bar],"black")
    hold on
    plot(xj-dx,int_L,"*red")
    hold on
    plot(xj,int_I,"*red")
    hold on
    plot(xj+dx,int_R,"*red")

    xlabel("x [arb].")
    ylabel("y [arb].")

    fprintf("\n");
end

end