function [x,f,iters,grad_norm] = newtmeth_eigmod(fun, x0, tol, maxit,c)

x(:,1) = x0;
[f(1),g,H] = fun(x0);
grad_norm(1) = norm(g);

% Print the initial value information
fprintf('Iter:%3d | function value: %4.5f \n',0,f(1));
fprintf('------------------------------------------------------------------------\n');
% Start the loop from 1 to maximum number of iterations. If the tolerance
% is reached this loop is terminated using the break command.
for iters = 1:maxit
    % Factorize the Hessian with modified chol
    l = x(1,iters);
    w = x(2,iters);
    lamda_min = 0.5*(l^2+w^2+2*c)-0.5*sqrt((l^2+w^2)^2+12*(w^2)*(l^2));
    if lamda_min>0
    	tau = 0;
    else
    	% tau = max(2*abs(lamda_min),1e-3);
    	tau = abs(lamda_min)+0.01;
    end
 %    tau = 0;
 %    while true;
	% [R,p] = chol(H+tau*eye(2,2));
	% if ~p
	% 	break;
	% else
	% 	tau = max([2*tau, 1e-3]);
	% end
	% end
    H(1,1) = H(1,1)+tau;
    H(2,2) = H(2,2)+tau;
    % Calculate the ascend direction as (H+tau*I)^{-1}*g with "\" operator
    % and the factorization result.
    dir = zeros(2,1);
    det_hess = H(1,1)*H(2,2)-H(1,2)*H(2,1);
    dir(1,1) = (1/det_hess)*(H(2,2)*g(1)-H(1,2)*g(2));
    dir(2,1) = (1/det_hess)*(-H(2,1)*g(1)+H(1,1)*g(2));
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
    % if f(iters+1)>f(iters)
    %     error('not descent step observed!');
    % end
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
