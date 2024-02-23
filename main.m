%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/20/2024
% IMEX Implementation Test:
% Verifying results of figure 8: 
% (https://arxiv.org/pdf/2102.11939.pdf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Setup the problem
app = make_app();

% Iterate over the time domain
while(app.grid.time < app.grid.t_max)

    % Push the app by one step
    app = update_app(app);
    
    % Call diagnostics
    app = diagnostics(app);

end