function [app] = make_app()
% This function sets up the grid structure
% and associated constants etc for the problem 


%%% START PROBLEM 1: Variable Knudsen Number %%%
% Periodic BC

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
% Velocity:
grid.v_max = 15;
grid.v_min = -grid.v_max;
grid.Nv = 150;
grid.v = linspace(grid.v_min,grid.v_max,grid.Nv);
% Time:
grid.t_min = 0.0;
grid.t_max = 0.5;
grid.dt = (1/24)*(grid.dx/grid.v_max);
grid.time = grid.t_min;
grid.NT = 1;
grid.time_vec = [grid.t_min];

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
grid.R = mod( linspace(1,Nx,Nx), Nx) + 1;
grid.L = mod( linspace(-1,Nx-2,Nx), Nx) + 1;

% Build phase-space grids:
app.f = zeros(grid.Nx,grid.Nv);
app.fp_im_half = zeros(grid.Nx,grid.Nv);
app.fm_ip_half = zeros(grid.Nx,grid.Nv);
app.fp_im_half_limit = zeros(grid.Nx,grid.Nv);
app.fm_ip_half_limit = zeros(grid.Nx,grid.Nv);
app.F_jp_half = zeros(grid.Nx,grid.Nv);
app.F_jm_half = zeros(grid.Nx,grid.Nv);

% Make the Knudsen Number array
eps0 = 10e-5;
grid.eps = problem1_Knudsen_num(eps0,grid.x);

% Save the grid object to the app 
app.grid = grid;

% Add any external app functions

%%% END PROBLEM 1: Vairable Knudsen Number %%%

end

% Knudsen Number
function eps = problem1_Knudsen_num(eps0,x)
    eps = eps0 + tanh(1 - 11*(x-1)) + tanh(1 + 11*(x-1)); 
end