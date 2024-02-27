%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/27/2024
% Test the Maxwellian Projection and correction routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Build required quantities:
grid.WENO_order = 3;

%%% WENO TESTING %%%
grid.moments_type = "Simple_No_Weno_reconst_fv";
grid.Nx = 1;
grid.Nv = 150;
Nx = grid.Nx;
Nv = grid.Nv;
f_IC = zeros(Nx,Nv);
f_fixed = zeros(Nx,Nv);
grid.x_max = 1;
grid.x_min = 0;
grid.Lx = grid.x_max - grid.x_min;
grid.x = linspace(grid.x_min,grid.x_max,grid.Nx);

grid.v_max = 15;
grid.v_min = -grid.v_max;
grid.Lv = grid.v_max - grid.v_min;
grid.v = linspace(grid.v_min,grid.v_max,grid.Nv);
grid.dv = grid.v(2) - grid.v(1);

% Save the grid to an app
app.grid = grid;

% Max testing
app.m0 = 1;
app.kb = 1;

% Intial moments:
n0 = 1;
u0 = 1; 
T0 = 1;

% Make a maxwellian distribtuion
for i = 1:grid.Nv
    f_IC(i) = maxwellian(n0,u0,T0,grid.v(i),app);
end

% Fake t for ref
t = [0,1];
ts = [0.25,0.75];

% Fix the maxwellian
f_fixed = fix_max(f_IC,n0,u0,T0,app);

% Initial projection/fixed
[n1,u1,T1] = moments(f_IC,app);
[n2,u2,T2] = moments(f_fixed,app);
n_corr = [n1,n2];
u_corr = [u1,u2];
T_corr = [T1,T2];


% Plot the exact solution
subplot(2,3,1)
plot(grid.v,f_IC,LineWidth=2)
hold on
plot(grid.v,f_fixed,":",LineWidth=2)
title("f(v)")
xlabel("v")
ylabel("f(v)")
legend("IC","Fixed")


% Plot the moments
subplot(2,3,2)
plot(t,n0 + t*0,"black")
hold on
plot(ts,n_corr,"*")
title("n")
xlabel("Left: IC, Right Fixed")
ylabel("n")
legend("Specified","Fixed")

% Plot the moments
subplot(2,3,3)
plot(t,u0 + t*0,"black")
hold on
plot(ts,u_corr,"*")
title("u")
xlabel("Left: IC, Right Fixed")
ylabel("u")
legend("Specified","Fixed")

% Plot the moments
subplot(2,3,4)
plot(t,T0 + t*0,"black")
hold on
plot(ts,T_corr,"*")
title("T")
xlabel("Left: IC, Right Fixed")
ylabel("T")
legend("Specified","Fixed")