function M_Eq = fix_max(M_Eq,n,u,T,app)
% fix Maxwellian distribtion

% Compute the moments 
[n_c,u_c,T_c] = moments(M_Eq,app);

% Simple fix of density only:
M_Eq = (n/n_c)*M_Eq;


% % Compute differences 
% n_diff = rel_diff(n_c,n);
% u_diff = rel_diff(u_c,u);
% T_diff = rel_diff(T_c,T);
% max_diff = max([n_diff,u_diff,T_diff]);
% 
% % Number of fixing iterations:
% niter = 0;
% tol = 1e-12;
% 
% % Iterate to correct the points
% while (niter < 15 && max_diff > tol)
% 
%     % Run the moment correction routine:
% 
%     % Compute the new distribution
% 
%     % Recompute moments 
% 
%     % Compute the error
%     
%     % Advance niter
%     niter = niter + 1;
% end

end




% Rel_diff
function res = rel_diff(a,b)

% Compute the relative difference between a and b
res = abs(a - b)/abs((a+b)/2);

end