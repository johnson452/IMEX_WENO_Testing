function [app] = update_app(app)
% Pushes the application by one timestep

% Grab the grid:
grid = app.grid;

% Time each push
tic

% Call the time-integration method:
if strcmp(grid.scheme,"JHU_FO") == 1
    app = SSP_IMEX_FIRST_ORDER(app);
elseif strcmp(grid.scheme,"JHU_SO") == 1
    app = SSP_IMEX_SECOND_ORDER(app);
else
    fprintf("No advancement method chosen!\n");
end

% End time
toc

% Push the time, save grid changes
grid.time = grid.time + grid.dt;
grid.NT = grid.NT + 1;
grid.time_vec = [grid.time_vec,grid.time];
app.grid = grid;

% Output the
fprintf("[%1.1f %%] Iteration (%d)\n",100*grid.time/grid.t_max,app.grid.NT)

end