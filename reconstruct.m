% Reconstruct
function [poly_array_tilde] = reconstruct(f,dir,app)

% WENO reconstruct -> gives f_{bar}, fp_{kp_half}, fm_{km_half}
fp_km_half = zeros(size(f));
fm_kp_half = zeros(size(f));
[fp_km_half,fm_kp_half,f_bar] = WENO(f,fp_km_half,fm_kp_half,dir,app.grid);
poly_array_tilde = WENO_poly_build_and_limit(fp_km_half,fm_kp_half,f_bar,dir,app.grid);

end
