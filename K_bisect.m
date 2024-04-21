function [soln, er_est] = K_bisect(R, xL, xR, tol, maxIter)

count = 1;

while count <= maxIter
    xG = (xL+xR)/2;
    er_est = abs(xR-xL)/2;
    fG = R(xG);
    if er_est > tol
        if min(fG,[],"all") >= 0 
            xR = xG;
        elseif min(fG,[],"all") < 0
            xL = xG;
        end
        count = count + 1;
    else
        break
    end
end

soln = xG;