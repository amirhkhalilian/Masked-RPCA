clc; clear; close all;

maxit = 20;
gam = 0.1;

X(:,1) = [0;0];

for iter = 1:maxit
Y(1,iter) = (2+X(1,iter)/gam)/(2+1/gam);
Y(2,iter) = X(2,iter);

Z(1,iter) = 2*Y(1,iter)-X(1,iter);
Z(2,iter) = (10*Z(1,iter)^2+(1/gam)*(2*Y(2,iter)-X(2,iter)))/(10+1/gam);

X(1,iter+1) = X(1,iter) + Z(1,iter)-Y(1,iter);
X(2,iter+1) = X(2,iter) + Z(2,iter)-Y(2,iter);

f(iter) = (Y(1,iter)-1)^2;
g(iter) = 5*(Z(2,iter)-Z(1,iter)^2)^2;
D(iter) = f(iter) + g(iter) - 1/(2*gam)*((Y(1,iter)-Z(1,iter))^2+(Y(2,iter)-Z(2,iter))^2) +...
 			(1/gam) *(((Z(1,iter)-Y(1,iter))*(X(1,iter)-Y(1,iter)))+((Z(2,iter)-Y(2,iter))*(X(2,iter)-Y(2,iter))));

fprintf('i:%2d|f:%1.4e|g:%1.4e|D:%1.4e\n',iter,f(iter),g(iter),D(iter));

end

