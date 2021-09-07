function [ alpha ] = BTLS( fun, x, p_k, c, rho )
% This is the Backtracking line search implemented 
% using the method in Alg3.1 N&W

% Inputs:
% fun: handle to function
% x: current point
% p_k: decent direction
% c: constant in (0,1)
% rho: constant in (0,1)

% Output:
% alpha: BTLS step-size


% check inputs
if (c>1 || c<0)
    error('Invalid input c');
end
if (rho>1 || rho<0)
    error('Invalid input rho');
end

% initialize alpha
alpha=1;

% while loop to find the step-size

while true
    % claculate the function and gradiant at x
    [f_x,g_x] = fun(x);
    [f_next] = fun(x+alpha*p_k);
    if f_next > f_x + c * alpha * g_x' * p_k
        alpha = rho * alpha;
    else
        break;
    end  
end
end

