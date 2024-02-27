% Map the the grid indices
function [indx_arr] = mapindex(dir,indx_arr,Nx)
if strcmp(dir,"x") == 1
    for i = 1:max(size(indx_arr))
        if (indx_arr(i) > Nx)
            indx_arr(i) = indx_arr(i) - Nx;
        end
        if (indx_arr(i) < 1)
            indx_arr(i) = indx_arr(i) + Nx;
        end
    end
elseif (strcmp(dir,"v") == 1)
    for i = 1:max(size(indx_arr))
        if (indx_arr(i) > Nx)
            indx_arr(i) = Nx;
        end
        if (indx_arr(i) < 1)
            indx_arr(i) = 1;
        end
    end
end
end