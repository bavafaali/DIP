function [kernel] = twoMatchedFilterKernel(k, sigma)
    L = 7;
    halfL = floor(L/2);

    r_i = [cos]    

    [x, y] = meshgrid(-halfL:halfL, -halfL:halfL);
    kernel = zeros(L, L);

    for i = 1:L
        for j = 1:L
            if abs(y(i,j)) <= L/2
                kernel(i,j) = -k * exp(-x(i,j)^2 / (2 * sigma^2));
            else
                kernel(i,j) = 0;
            end
        end
    end

    kernel = kernel / sum(abs(kernel(:)));
    disp(kernel);
end
