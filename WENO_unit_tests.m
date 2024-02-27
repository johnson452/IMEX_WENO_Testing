%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/21/2024
% Testing for the code
%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%Order for testing
for k = 1:3
    grid.WENO_order = k;
    grid.moments_type = "WENO_Reconstructed_fv";


    %%% WENO TESTING %%%
    grid.Nx = 7;
    grid.Nv = 1;
    Nx = grid.Nx;
    Nv = grid.Nv;
    f_full = zeros(Nx,Nv);
    fp_im_half = zeros(Nx,Nv);
    fm_ip_half = zeros(Nx,Nv);
    fp_im_half_limit = zeros(Nx,Nv);
    fm_ip_half_limit = zeros(Nx,Nv);
    grid.x_max = 1;
    grid.x_min = 0;
    grid.Lx = grid.x_max - grid.x_min;
    grid.x = linspace(grid.x_min,grid.x_max,grid.Nx);
    grid.dx = grid.x(2) - grid.x(1);
    %grid.R = mod( linspace(1,Nx,Nx), Nx) + 1;
    %grid.L = mod( linspace(-1,Nx-2,Nx), Nx) + 1;
    %grid.I = linspace(1,Nx,Nx);
    %grid.xl = grid.x - grid.dx/2;
    %grid.xr = grid.x + grid.dx/2;

    % Bumpy test:
    % make inital f_full
    % for i = 1:Nx
    %     if i < 10
    %         f_full(i) = 0;
    %     elseif (i >= 10) && (i <= 20)
    %         f_full(i) = 1 + 2*sin(8*pi*grid.x(i)/grid.Lx);
    %     elseif i > 20
    %         f_full(i) = -exp(-grid.x(i)/grid.Lx);
    %     end
    % end

    % Smooth test:
    for i = 1:Nx
        f_full(i) = 1.0 + sin(2*pi*grid.x(i)/(grid.Lx+grid.dx));
    end

    %Discontinuity:
    % for i = 1:Nx
    %     if grid.x(i) > 0.5
    %         f_full(i) = 1;
    %     else
    %         f_full(i) = 0;
    %     end
    % end

% 
%     % Velocity Smooth test:
%     for i = 1:Nv
%         f_full_v(i) = sin(2*pi*grid.v(i)/(grid.Lv+grid.dv));
%     end


    % Plot the exact solution:
    Nx_exact = 10000;
    f_exact = zeros(1,Nx_exact);
    x_exact = linspace(grid.x_min,grid.x_max,Nx_exact);
    dx_exact = x_exact(2) - x_exact(1);
    for i = 1:Nx_exact
        f_exact(i) = 1 + sin(2*pi*x_exact(i)/(grid.Lx+grid.dx));
    end


    % Call WENO
    dir = "x";
    [fp_im_half,fm_ip_half,f_bar] = WENO(f_full,fp_im_half,fm_ip_half,dir,grid);
    [fp_im_half_limit,fm_ip_half_limit] = WENO_flux_limit(fp_im_half,fm_ip_half,f_bar,fp_im_half_limit,fm_ip_half_limit,grid);
    poly_array_tilde = WENO_poly_build_and_limit(fp_im_half,fm_ip_half,f_bar,dir,grid);


    % Plot
    subplot(4,3,k)
    plot(x_exact,f_exact,"color","black")
    hold on
    plot(grid.x,f_full,"*-")
    hold on
    plot(grid.x,f_bar,"*")
    legend("f_exact","f_{full}","f_{bar}")
    title(sprintf("f_{bar}. (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("x")
    ylabel("f(x)")

    subplot(4,3,k+3)
    plot(x_exact,f_exact,"color","black")
    hold on
    plot(grid.x,f_full,"*-")
    hold on
    plot(grid.x - grid.dx/2,fp_im_half,"*")
    hold on
    plot(grid.x + grid.dx/2,fm_ip_half,"*")
    title(sprintf("WENO Interp. (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("x")
    ylabel("f(x)")
    legend("f_exact","f_{full}","f_{i-1/2}^{+}","f_{i+1/2}^{-}")

    subplot(4,3,k+6)
    plot(x_exact,f_exact,"color","black")
    hold on
    plot(grid.x,f_full,"*-")
    hold on
    plot(grid.x - grid.dx/2,fp_im_half_limit,"*")
    hold on
    plot(grid.x + grid.dx/2,fm_ip_half_limit,"*")
    title(sprintf("WENO Interp. Limited (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("x")
    ylabel("f(x)")
    legend("f_exact","f_{full}","f_{i-1/2}^{+}","f_{i+1/2}^{-}")

    subplot(4,3,k+9)
    plot(x_exact,f_exact,"color","black")
    hold on
    plot(grid.x,f_full,"*-")
    hold on
    for i = 1:Nx
        local_x_vec = linspace(grid.x(i)-grid.dx/2,grid.x(i)+grid.dx/2,20);
        local_poly_evaled = poly_eval(poly_array_tilde(:,i,1),grid.x(i),local_x_vec);
        plot(local_x_vec,local_poly_evaled,"--")
    end
    title(sprintf("WENO Interp. Limited (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("x")
    ylabel("f(x)")

end



%Order for testing
figure2 = figure();
for k = 1:3
    grid.WENO_order = k;

    %%% WENO TESTING %%%
    grid.Nx = 3;
    grid.Nv = 7;
    Nx = grid.Nx;
    Nv = grid.Nv;
    f_full = zeros(Nx,Nv);
    fp_im_half = zeros(Nx,Nv);
    fm_ip_half = zeros(Nx,Nv);
    fp_im_half_limit = zeros(Nx,Nv);
    fm_ip_half_limit = zeros(Nx,Nv);
    grid.x_max = 1;
    grid.x_min = 0;
    grid.Lx = grid.x_max - grid.x_min;
    grid.x = linspace(grid.x_min,grid.x_max,grid.Nx);
    %grid.dx = grid.x(2) - grid.x(1);

    grid.v_max = 1;
    grid.v_min = -1;
    grid.Lv = grid.v_max - grid.v_min;
    grid.v = linspace(grid.v_min,grid.v_max,grid.Nv);
    grid.dv = grid.v(2) - grid.v(1);
    %grid.R = mod( linspace(1,Nx,Nx), Nx) + 1;
    %grid.L = mod( linspace(-1,Nx-2,Nx), Nx) + 1;
    %grid.I = linspace(1,Nx,Nx);
    %grid.xl = grid.x - grid.dx/2;
    %grid.xr = grid.x + grid.dx/2;

    % Bumpy test:
    % make inital f_full
    % for i = 1:Nx
    %     if i < 10
    %         f_full(i) = 0;
    %     elseif (i >= 10) && (i <= 20)
    %         f_full(i) = 1 + 2*sin(8*pi*grid.x(i)/grid.Lx);
    %     elseif i > 20
    %         f_full(i) = -exp(-grid.x(i)/grid.Lx);
    %     end
    % end

    % Smooth test:
    %for i = 1:Nx
    %    f_full(i) = sin(2*pi*grid.x(i)/(grid.Lx+grid.dx));
    %end

    %Discontinuity:
    % for i = 1:Nx
    %     if grid.x(i) > 0.5
    %         f_full(i) = 1;
    %     else
    %         f_full(i) = 0;
    %     end
    % end


    % Velocity Smooth test:
    for j = 1:Nx
    for i = 1:Nv
        f_full(j,i) = sin(2*pi*grid.v(i)/(grid.Lv+grid.dv));
    end
    end


    % Plot the exact solution:
    Nv_exact = 10000;
    f_exact = zeros(1,Nv_exact);
    v_exact = linspace(grid.v_min,grid.v_max,Nv_exact);
    dv_exact = v_exact(2) - v_exact(1);
    for i = 1:Nv_exact
        f_exact(i) = sin(2*pi*v_exact(i)/(grid.Lv+grid.dv));
    end


    % Call WENO
    dir = "v";
    [fp_im_half,fm_ip_half,f_bar] = WENO(f_full,fp_im_half,fm_ip_half,dir,grid);
    poly_array_tilde = WENO_poly_build_and_limit(fp_im_half,fm_ip_half,f_bar,dir,grid);

    % Grab slice:
    f_full = f_full(2,:);
    poly_array_tilde = poly_array_tilde(:,2,:);
    fp_im_half = fp_im_half(2,:);
    fm_ip_half = fm_ip_half(2,:);
    f_bar = f_bar(2,:);

    % Plot
    subplot(3,3,k)
    plot(v_exact,f_exact,"color","black")
    hold on
    plot(grid.v,f_full,"*-")
    hold on
    plot(grid.v,f_bar,"*")
    legend("f_exact","f_{full}","f_{bar}")
    title(sprintf("f_{bar}. (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("x")
    ylabel("f(x)")

    subplot(3,3,k+3)
    plot(v_exact,f_exact,"color","black")
    hold on
    plot(grid.v,f_full,"*-")
    hold on
    plot(grid.v - grid.dv/2,fp_im_half,"*")
    hold on
    plot(grid.v + grid.dv/2,fm_ip_half,"*")
    title(sprintf("WENO Interp. (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("v")
    ylabel("f(v)")
    legend("f_exact","f_{full}","f_{i-1/2}^{+}","f_{i+1/2}^{-}")

    subplot(3,3,k+6)
    plot(v_exact,f_exact,"color","black")
    hold on
    plot(grid.v,f_full,"*-")
    hold on
    for i = 1:Nv
        local_v_vec = linspace(grid.v(i)-grid.dv/2,grid.v(i)+grid.dv/2,20);
        local_poly_evaled = poly_eval(poly_array_tilde(:,1,i),grid.v(i),local_v_vec);
        plot(local_v_vec,local_poly_evaled,"--")
    end
    title(sprintf("WENO Interp. Limited (Order: %d) (Poly Order: %d)",2*k -1,k-1))
    xlabel("v")
    ylabel("f(v)")
    %legend("f_exact","f_{full}","f_poly_rep")

end