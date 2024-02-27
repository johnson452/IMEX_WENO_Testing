%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 2/27/2024
% Check fourth order polyfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Allocate the syms
syms a0 a1 a2 a3 a4  dx x xj   uj ujm ujp u_jp_half u_jm_half  nu

% Create the poly function
poly = a0 + a1*(x-xj) + a2*(x-xj)^2 + a3*(x-xj)^3 + a4*(x-xj)^4;

% Create the constraint equations
eqn1_LHS = simplify(subs(poly, x,xj - dx/2));  % = u_jm_half
eqn2_LHS = simplify(subs(poly, x,xj + dx/2));  % = u_jp_half

% Integrate the polynomial
int_poly = int(poly);
eqn3_LHS = simplify((1/dx) * ( subs(int_poly,xj + 3*dx/2) - subs(int_poly,xj + dx/2) )); % = ujp
eqn4_LHS = simplify((1/dx) * ( subs(int_poly,xj + dx/2) - subs(int_poly,xj - dx/2) )); % = uj
eqn5_LHS = simplify((1/dx) * ( subs(int_poly,xj - dx/2) - subs(int_poly,xj - 3*dx/2) )); % = ujm


% Make the constraint equations
M = sym(zeros(5,5)); 
a = [a0,a1,a2,a3,a4];
b = [u_jm_half;u_jp_half;ujp;uj;ujm];
M(1,:) = coeffs(eqn1_LHS,a);
M(2,:) = coeffs(eqn2_LHS,a);
M(3,:) = coeffs(eqn3_LHS,a);
M(4,:) = subs(coeffs(eqn4_LHS + nu*a3 + nu*a1,a),nu,0);
M(5,:) = coeffs(eqn5_LHS,a);
aT = [a4;a3;a2;a1;a0];
Ma = M*aT;
 
% Verify Ma matches:
% disp(Ma(1) - eqn1_LHS);
% disp(Ma(2) - eqn2_LHS);
% disp(Ma(3) - eqn3_LHS);
% disp(Ma(4) - eqn4_LHS);
% disp(Ma(5) - eqn5_LHS);

% Compute M_inv*b
a_solve = simplify(inv(M)*b);