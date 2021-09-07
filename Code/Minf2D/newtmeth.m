function [x,f,iters,grad_norm] = newtmeth(fun, x0, tol, maxit)
% This function implements the newton method using backtracking
% line search and Cholesky with added Multiple of the identity. 
%
% Inputs:
% fun: function handle
% x0: initial point
% tol: tolerance on lambda^2
% maxit: maximum number of iterations
%
% Outputs:
% x: points at each step
% f: function values at each step
% iters: number of iterations taken
% grad_norm: norm of the gradient at each iteration

% Compute the function and gradient and Hessian at initial point
[f(1),g,H] = fun(x0);
x(:,1) = x0;    % matrix x will save all the x vectors calculated
grad_norm(1) = norm(g);
% Print the initial value information
fprintf('Iter:%3d | function value: %4.5f \n',0,f(1));
fprintf('------------------------------------------------------------------------\n');
% Start the loop from 1 to maximum number of iterations. If the tolerance
% is reached this loop is terminated using the break command.
for iters = 1:maxit
    % Factorize the Hessian with modified chol
    [R,tau] = chol_add_mult_I(H,1e-3);
    % Calculate the ascend direction as (H+tau*I)^{-1}*g with "\" operator
    % and the factorization result.
    dir = R\(R'\g);
    % Calculate t using BTLS. direction is -dir
    t = BTLS(fun,x(:,iters),-dir,0.25,0.5);
    % Calculate the next step x vector and save it in matrix x.
    x(:,iters+1) = x(:,iters)+t*(-dir);
    % Calculate the function and gradient an Hessian at next point
    [f(iters+1),g,H] = fun(x(:,iters+1));
    % Calculate the norm of the gradient to plot later
    grad_norm(iters+1) = norm(g);
    % Since newton method is descend algorithm if the value of function is
    % increased an error is generated.
    if f(iters+1)>f(iters)
        error('not descent step observed!');
    end
    % print the iteration results
    fprintf('Iter: %4d | t: %2.4f | grad_norm: %1.2e| function value: %2.4f| tau:%1.2e \n',iters,t,grad_norm(iters+1),f(iters+1),tau);
    % Check the stopping criterion
    if norm(g) < tol
        break
    end
end
% Print summary of the results.
fprintf('------------------------------------------------------------------------\n');
fprintf('Summary of Results\n');
fprintf('------------------------------------------------------------------------\n');
if iters<maxit
    fprintf('tolerance reached.\n');
    fprintf('norm of the gradient = %1.4e\n',norm(g));
    fprintf('Optimal value = %2.4f\n',f(end));
elseif iters==maxit
    fprintf('maximum iterations reached without reaching the tolerance\n');
    fprintf('norm of the gradient = %1.4e\n',norm(g));
    fprintf('final function value = %2.4f\n',f(end));
end
fprintf('------------------------------------------------------------------------\n');

end

