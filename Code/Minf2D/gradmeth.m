function [x,f,norm_grad,iters] = gradmeth(fun, x0, tol, maxit)
% This function implements the gradient decent method using backtracking
% line search. 
%
% Inputs:
% fun: function handle
% x0: initial point
% tol: tolerance on norm of the gradient
% maxit: maximum number of iterations
%
% Outputs:
% x: points at each step
% f: function values at each step
% iters: number of iterations taken

% Compute the function and gradient and Hessian at initial point
[f(1),g] = fun(x0);
x(:,1) = x0;        % matrix x will save all the x vectors calculated 
norm_grad(1) = norm(g);
% Display the initial point information
fprintf('Iter:%3d | function value: %4.5f \n',0,f(1));
fprintf('--------------------------------------------------------------\n'); 
% Start the loop from 1 to maximum number of iterations. If the tolerance
% is reached this loop is terminated using the break command.
for iters = 1:maxit
    % find the step-size using BTLS with c = 1e-4 rho = 0.5.
    t = BTLS(fun,x(:,iters),-g,1e-4,0.5);
    % change the x in the -g direction with the step-size t
    x(:,iters+1) = x(:,iters)+t*(-g);
    % calculate the function and the gradient at next point
    [f(iters+1),g] = fun(x(:,iters+1));
    norm_grad(iters+1) = norm(g);
    % print the step information
    % fprintf('Iter: %4d | alpha: %1.2e | f_val: %2.4f | |g|: %1.5e \n',iters,t,f(iters+1),norm_grad(iters+1));
    % check convergance
    if norm_grad(iters+1)<tol
        break
    end
    % since Gradient descent always decreases the function an error is
    % generated if increase is observed
    if f(iters+1)>f(iters)
        error('not descent step observed!');
    end
end
% summary of the results 
fprintf('--------------------------------------------------------------\n');
fprintf('Summary of Results\n');
fprintf('--------------------------------------------------------------\n');
if iters<maxit
    fprintf('tolerance reached.\n');
    fprintf('norm of the gradient = %1.4e\n',norm(g));
    fprintf('Optimal value = %2.4f\n',f(end));
elseif iters==maxit
    fprintf('maximum iterations reached without reaching the tolerance\n');
    fprintf('norm of the gradient = %1.4e\n',norm(g));
    fprintf('final function value = %2.4f\n',f(end));
end
fprintf('--------------------------------------------------------------\n');    

end
