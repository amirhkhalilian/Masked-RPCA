function [x,ng] = stupid_solver(x0,x_l,x_w,tol, maxit, c);
	x(:,1) = x0;
	g1 = (x(2,1))^2*(x(1,1)) + c*((x(1,1))-x_l);
	g2 = (x(1,1))^2*(x(2,1)) + c*((x(2,1))-x_w);
	ng(1) = norm([g1;g2]);
	for iter=1:maxit
		x(1,iter+1) = (c*x_l)./(x(2,iter)^2+c);
		x(2,iter+1) = (c*x_w)./(x(1,iter+1)^2+c);
		g1 = (x(2,iter+1))^2*(x(1,iter+1)) + c*((x(1,iter+1))-x_l);
		g2 = (x(1,iter+1))^2*(x(2,iter+1)) + c*((x(2,iter+1))-x_w);
		ng(iter+1) = norm([g1;g2]);
		fval =0;
		fprintf('i:%2d|ng:%1.2e|f:%1.2e\n',iter,ng(iter+1),fval);
		if ng(iter+1)<=tol
			break;
		end
	end
end