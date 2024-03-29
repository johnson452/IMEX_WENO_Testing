function [app] = make_app()
% This function sets up the grid structure
% and associated constants etc for the problem 


%%% START PROBLEM 1: Variable Knudsen Number %%%
% Periodic BC
grid.WENO_order = 3;
grid.moments_type = "Simple_No_Weno_reconst_fv"; %,"WENO_Reconstructed_fv";
grid.scheme = "JHU_FO"; %"SSP_RK3_EXPLICIT"; % "JHU_SO"; % "JHU_FO"

% Constants
app.m0 = 1;
app.kb = 1;

% Build the grid object
% Spatial:
grid.Nx = 40;
grid.x_min = 0;
grid.x_max = 2;
grid.x = linspace(grid.x_min,grid.x_max,grid.Nx);
grid.dx = grid.x(2) - grid.x(1);
grid.Lx = grid.x_max - grid.x_min;
% Velocity:
grid.v_max = 15;
grid.v_min = -grid.v_max;
grid.Nv = 150;
grid.v = linspace(grid.v_min,grid.v_max,grid.Nv);
grid.dv = grid.v(2) - grid.v(1);
% Time:
grid.t_min = 0.0;
grid.t_max = 0.5;
grid.dt = (1/48)*(grid.dx/grid.v_max); %(1/24)*(grid.dx/grid.v_max); %(1/240)*(grid.dx/grid.v_max); 
grid.time = grid.t_min;
grid.NT = 1;
grid.time_vec = [grid.t_min];
grid.diag_interval  = 30;
fprintf("We expect the simulation to take: %d iterations\n",ceil(grid.t_max/grid.dt));


% Build velocity grids:
grid.v_plus = zeros(grid.Nv,1);
grid.v_minus = zeros(grid.Nv,1);
for i = 1:grid.Nv
    if grid.v(i) > 0
        grid.v_plus(i) = grid.v(i);
        grid.v_minus(i) = 0;
    else
        grid.v_plus(i) = 0;
        grid.v_minus(i) = grid.v(i);
    end
end

% Build utility arrays
grid.R = mod( linspace(1,grid.Nx,grid.Nx), grid.Nx) + 1;
grid.L = mod( linspace(-1,grid.Nx-2,grid.Nx), grid.Nx) + 1;

% Build phase-space grids:
app.f = zeros(grid.Nx,grid.Nv);
app.fp_im_half = zeros(grid.Nx,grid.Nv);
app.fm_ip_half = zeros(grid.Nx,grid.Nv);
app.fp_im_half_limit = zeros(grid.Nx,grid.Nv);
app.fm_ip_half_limit = zeros(grid.Nx,grid.Nv);
app.F_jp_half = zeros(grid.Nx,grid.Nv);
app.F_jm_half = zeros(grid.Nx,grid.Nv);

% % Make IC - Knuudsen Number Problem
for i = 1:grid.Nx
    for j = 1:grid.Nv
        
        % Grid location:
        x = grid.x(i);
        v = grid.v(j);

        % Moments:
        sin_term = 0.2*sin(pi*x*((grid.Lx-grid.dx)/grid.Lx));
        n = 1 + sin_term;
        u = 1;
        T = 1/(1 + sin_term);
        
        % Biuld f:
        %app.f(i,j) = maxwellian(n,u,T,v,app);
        app.f(i,j) = 0.5*maxwellian(n,u,T,v,app) + 0.3*maxwellian(n,-0.5*u,T,v,app);
    end
end

% % Make IC - Bi-Maxwellian
% for i = 1:grid.Nx
%     for j = 1:grid.Nv
%         
%         % Grid location:
%         x = grid.x(i);
%         v = grid.v(j);
% 
%         % Moments:
%         n = 1;
%         u = 1;
%         T = 1;
%         
%         % Biuld f:
%         app.f(i,j) = maxwellian(n,u,T,v,app) + maxwellian(n,-u,T,v,app);
%         %app.f(i,j) = 0.5*maxwellian(n,u,T,v,app) + 0.3*maxwellian(n,-0.5*u,T,v,app);
%     end
% end

% Make IC - Phase Mixing
% for i = 1:grid.Nx
%     for j = 1:grid.Nv
%         
%         % Grid location:
%         x = grid.x(i);
%         v = grid.v(j);
% 
%         % Moments:
%         n = 1;
%         u = 1;
%         T = 1;
%         
%         % Biuld f:
%         sin_term = 0.2*sin(4*pi*x*((grid.Lx-grid.dx)/grid.Lx));
%         app.f(i,j) = 0.001*maxwellian(n,u,T,v,app);
%         
%         % Add inital distribution
%         if abs(v) < 7 && sin_term > 0
%             app.f(i,j) = app.f(i,j) + 1;
%         end
%     end
% end


% Make the Knudsen Number array
app.eps0 = 1e-5;

% Save the grid object to the app 
app.grid = grid;

% Save the inital distribution:
app.f_IC = app.f;
[n,u,T] = moments(app.f_IC,app);
app.n_IC = n;
app.u_IC = u;
app.T_IC = T;


% Add any external app functions

%%% END PROBLEM 1: Vairable Knudsen Number %%%

end