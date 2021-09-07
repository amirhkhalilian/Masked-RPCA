clc; clear; close all;

f = @(w,l) (w^2)*(l^2);
g = @(w,l) [2*l*w^2;2*w*l^2];

x = [2;2];rand(2,1);%[1;1];
u = [0;1];
counter = 1;
for h = 10.^[-15:1:0]
	ratio(counter) = norm(g(x(1)+h*u(1),x(2)+h*u(2))-g(x(1),x(2)))/norm(h*u);
	counter = counter+1;
end
2*sqrt((8*x(1)^2*x(2)^2+x(1)^4+x(2)^4+4*x(1)*x(2)^3+4*x(2)*x(1)^3)/2)
ratio(5)

loglog(10.^[-15:1:0],ratio);