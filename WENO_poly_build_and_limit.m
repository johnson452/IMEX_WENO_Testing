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
x = grid.x;

% If the direction is v:
if strcmp(dir,"v") == 1
    Nx = sz_f(2);
    dx = grid.dv;
    x = grid.v;
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
        jm = mapindex(dir,j-1,Nx);
        jp = mapindex(dir,j+1,Nx);

        % Compute poly coeff
        a0 = (1/192)*(f_bar(jm,v_index) + 298*f_bar(j,v_index) + f_bar(jp,v_index) - 54*(fp_im_half(j,v_index) + fm_ip_half(j,v_index)));
        a1 = (1/(8*dx))*(f_bar(jm,v_index) - f_bar(jp,v_index) - 10*(fp_im_half(j,v_index) - fm_ip_half(j,v_index)));
        a2 = (1/(8*dx^2))*( -(f_bar(jm,v_index) + 58*f_bar(j,v_index) + f_bar(jp,v_index)) + 30*(fp_im_half(j,v_index) + fm_ip_half(j,v_index)));
        a3 = (1/(2*dx^3))*(f_bar(jp,v_index) - f_bar(jm,v_index) +2*(fp_im_half(j,v_index) - fm_ip_half(j,v_index)));
        a4 = (1/(12*dx^4))*( (5*f_bar(jm,v_index) + 50*f_bar(j,v_index) + 5*f_bar(jp,v_index)) - 30*(fp_im_half(j,v_index) + fm_ip_half(j,v_index)));

        % Poly eval:
        poly = [a0,a1,a2,a3,a4];

        % Check polynomial matches conditions of:
        % (Sec 2.4: https://www.sciencedirect.com/science/article/pii/S0021999109007165)
        %check_poly_fourth_order(poly,f_bar(j,v_index),f_bar(jm,v_index),...
        %    f_bar(jp,v_index),fp_im_half(j,v_index),fm_ip_half(j,v_index),x(j),dx);

        % (3.35a,b) of: https://epubs.siam.org/doi/epdf/10.1137/17M1144362
        xi_j = (f_bar(j,v_index) - fp_im_half(j,v_index)*w1 - fm_ip_half(j,v_index)*w4)/(w2 + w3);
        xj = x(j);
        mj = min([poly_eval(poly,xj,xj-dx/2),poly_eval(poly,xj,xj+dx/2),xi_j]); % equiv to eval_poly, see PP before
        sigma_j = sigma(mj,f_bar(j,v_index));

        %Limit the polynomial
        poly_tilde = sigma_j*(poly) - sigma_j*[1,0,0,0,0]*f_bar(j,v_index) + [1,0,0,0,0]*f_bar(j,v_index);

        % Save the polynomial
        poly_array_tilde(:,j,v_index) = poly_tilde;

        % Check the limited version
        %check_poly_fourth_order_limited(poly_tilde,f_bar(j,v_index),f_bar(jm,v_index),...
        %    f_bar(jp,v_index),fp_im_half(j,v_index),fm_ip_half(j,v_index),x(j),dx);

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



% Check the polynomial fourth order matches the condtions
function check_poly_fourth_order(poly,f_bar,fL_bar,fR_bar,fp_im_half,fm_ip_half,xj,dx)

% Does the edge evaluations match?
x_jm_half = xj - dx/2;
x_jp_half = xj + dx/2;
pj_x_jm_half = poly_eval(poly,xj,x_jm_half);
pj_x_jp_half = poly_eval(poly,xj,x_jp_half);
err_pj_x_jm_half = rel_diff(pj_x_jm_half,fp_im_half);
err_pj_x_jp_half = rel_diff(pj_x_jp_half,fm_ip_half);

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


if err_int > tol || err_pj_x_jm_half > tol || err_pj_x_jp_half > tol

    fprintf("(POLY EDGES) poly_comp(-): %1.16e, val(-): %1.16e, rel_diff(-): %1.16e\n",pj_x_jm_half,fp_im_half,err_pj_x_jm_half);
    fprintf("(POLY EDGES) poly_comp(+): %1.16e, val(+): %1.16e, rel_diff(+): %1.16e\n",pj_x_jp_half,fm_ip_half,err_pj_x_jp_half);
    fprintf("(POLY MEAN) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_I,f_bar,err_int);
    fprintf("(POLY MEAN L) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_L,fL_bar,err_L_int);
    fprintf("(POLY MEAN R) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_R,fR_bar,err_R_int);
    fprintf("\n");

    clf()
    x_lower = xj - 3*dx/2;
    x_upper = xj + 3*dx/2;
    x_span = linspace(x_lower,x_upper);
    y_span = poly_eval(poly,xj,x_span);
    plot(x_span,y_span,"red")
    hold on
    plot([xj - dx/2,xj + dx/2],[fp_im_half,fm_ip_half],"*black")
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


% Check the polynomial fourth order matches the condtions
function check_poly_fourth_order_limited(poly,f_bar,fL_bar,fR_bar,fp_im_half,fm_ip_half,xj,dx)

% Does the edge evaluations match?
x_jm_half = xj - dx/2;
x_jp_half = xj + dx/2;
pj_x_jm_half = poly_eval(poly,xj,x_jm_half);
pj_x_jp_half = poly_eval(poly,xj,x_jp_half);
err_pj_x_jm_half = rel_diff(pj_x_jm_half,fp_im_half);
err_pj_x_jp_half = rel_diff(pj_x_jp_half,fm_ip_half);

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


% Gauss-Labatto Quad wieghts
w1 = 1/12;  % -1/2
w2 = 5/12;  % -0.447214/2
w3 = 5/12;  % +0.447214/2
w4 = 1/12;  % +1/2
xi_tilde = (f_bar - pj_x_jm_half*w1 - pj_x_jp_half*w4)/(w2 + w3);

tol = 1e-8;

if err_int > tol || pj_x_jm_half < 0 || pj_x_jp_half < 0 || xi_tilde < 0

    fprintf("(POLY EDGES) poly_comp(-): %1.16e, val(-): %1.16e, rel_diff(-): %1.16e\n",pj_x_jm_half,fp_im_half,err_pj_x_jm_half);
    fprintf("(POLY EDGES) poly_comp(+): %1.16e, val(+): %1.16e, rel_diff(+): %1.16e\n",pj_x_jp_half,fm_ip_half,err_pj_x_jp_half);
    fprintf("(POLY MEAN) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_I,f_bar,err_int);
    fprintf("(POLY MEAN L) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_L,fL_bar,err_L_int);
    fprintf("(POLY MEAN R) poly_comp: %1.16e, val: %1.16e, rel_diff: %1.16e\n",int_R,fR_bar,err_R_int);
    fprintf("\n");

    clf()
    x_lower = xj - 3*dx/2;
    x_upper = xj + 3*dx/2;
    x_span = linspace(x_lower,x_upper);
    y_span = poly_eval(poly,xj,x_span);
    plot(x_span,y_span,"red")
    hold on
    plot([xj - dx/2,xj + dx/2],[fp_im_half,fm_ip_half],"*black")
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