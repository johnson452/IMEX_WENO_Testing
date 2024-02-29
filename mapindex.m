% Map the the grid indices
function [indx_arr] = mapindex(dir,indx_arr,Nx)
sz = max(size(indx_arr));
    for i = 1:sz
        if (indx_arr(i) > Nx)
            indx_arr(i) = indx_arr(i) - Nx;
        end
        if (indx_arr(i) < 1)
            indx_arr(i) = indx_arr(i) + Nx;
        end
    end
end