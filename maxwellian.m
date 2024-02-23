function [val] =  maxwellian(n,u,T,v,app)
% Creates the value of the maxwellian at the point in space
% 1D only

% Grab constants
m = app.m0;
kb = app.kb;

% Create the value
val = n*sqrt(m/(2*pi*kb*T))*exp(-(m*(v-u)^2)/(2*kb*T));

end