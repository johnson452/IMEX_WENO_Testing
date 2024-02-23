function [y] = poly_eval(poly,xj,Sj)
% %Poly Eval:

% Size of solutions
sz = max(size(Sj));
sz_poly = max(size(poly));
y = zeros(sz,1);
for i = 1:sz
    for j = 1:sz_poly
        y(i) = y(i) + poly(j)*(Sj(i)-xj)^(j-1);
    end
end

end