%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/27/2024
% Check fourth order polyfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Allocate the syms
syms a0 a1 a2  dx x xj   uj ujm ujp  nu

% Create the poly function
poly = a0 + a1*(x-xj) + a2*(x-xj)^2;

% Create the constraint equations
% Integrate the polynomial
int_poly = int(poly);
eqn1_LHS = simplify((1/dx) * ( subs(int_poly,xj + 3*dx/2) - subs(int_poly,xj + dx/2) )); % = ujp
eqn2_LHS = simplify((1/dx) * ( subs(int_poly,xj + dx/2) - subs(int_poly,xj - dx/2) )); % = uj
eqn3_LHS = simplify((1/dx) * ( subs(int_poly,xj - dx/2) - subs(int_poly,xj - 3*dx/2) )); % = ujm


% Make the constraint equations
M = sym(zeros(3,3)); 
a = [a0,a1,a2];
b = [ujp;uj;ujm];
M(1,:) = coeffs(eqn1_LHS,a);
M(2,:) = subs(coeffs(eqn2_LHS + nu*a0 + nu*a1 + nu*a2,a),0); %subs(coeffs(eqn2_LHS + nu*a3 + nu*a1,a),nu,0);
M(3,:) = coeffs(eqn3_LHS,a);
aT = [a2;a1;a0];
Ma = M*aT;
 
% Verify Ma matches:
disp(Ma(1) - eqn1_LHS);
disp(Ma(2) - eqn2_LHS);
disp(Ma(3) - eqn3_LHS);

% Compute M_inv*b
a_solve = simplify(inv(M)*b);