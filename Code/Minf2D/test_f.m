clc; clear; close all;

x_l = 0.1;
x_w = 0.1;
gam = 2;
rho = 2;
maxL = 4;
maxW = 4;

% x_l = 1;
% x_w = 1;
% gam = 2;
% rho = 2;
% maxL = 4;
% maxW = 4;

% x_l = 3;
% x_w = 10;
% gam = 2;
% rho = 2;
% maxL = 15;
% maxW = 15;


minC = -0.5*(maxL^2+maxW^2)+sqrt(3*(maxL^2)*(maxW^2)+0.25*(maxL^2+maxW^2)^2);

fprintf('minimmum C value is: %+1.3f\n chosen c is: %+1.3f\n',minC,1/(rho*gam));

fun = @(x)f_eval_wl(x,x_l,x_w,rho,gam);

[z_nmeth_1,f_nmeth,iters,norm_grad_1] = newtmeth(fun, [2;-4], 1e-9, 100);
[z_nmeth_2,f_nmeth,iters,norm_grad_2] = newtmeth_eigmod(fun, [2;-4], 1e-9, 100,1/(gam*rho));
[z_nmeth_3,f_nmeth,iters,norm_grad_3] = newtmeth_eigabs(fun, [2;-4], 1e-9, 100);
[z_nmeth_4,norm_grad_4] = stupid_solver([2;-4],x_l,x_w,1e-9, 100,1/(gam*rho));
% [f,g,H] = f_eval_wl(z_nmeth_1(:,end),x_l,x_w,rho,gam);
% eig(H)

plot_2D_contour(maxL,maxW,x_l,x_w,rho,gam)
hold on;
plot(z_nmeth_1(1,:),z_nmeth_1(2,:),'r-*','MarkerSize',10);
plot(z_nmeth_2(1,:),z_nmeth_2(2,:),'b-*','MarkerSize',10);
plot(z_nmeth_3(1,:),z_nmeth_3(2,:),'m-*','MarkerSize',10);
plot(z_nmeth_4(1,:),z_nmeth_4(2,:),'k-*','MarkerSize',10);
legend('level sets','diag loading','add min eig','abs eig','Gauss-Seidel');
set(gca,'FontSize',14);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

figure;
semilogy(norm_grad_1,'r-o','LineWidth',2,'MarkerSize',10); hold on;
semilogy(norm_grad_2,'b-d','LineWidth',2,'MarkerSize',10);
semilogy(norm_grad_3,'m-s','LineWidth',2,'MarkerSize',10);
semilogy(norm_grad_4,'k-*','LineWidth',2,'MarkerSize',10);
grid on;
legend('diag loading','add min eig','abs eig','Gauss-Seidel');
xlabel('iteration number'); ylabel('gradient norm');
set(gca,'FontSize',14);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
% plot(z_nmeth(1,:),z_nmeth(2,:),'-*')

function plot_2D_contour(maxL,maxW,x_l,x_w,rho,gamma)
c = 1/(rho*gamma);
l = -maxL:0.01:maxL;
w = -maxW:0.01:maxW;
[L,W] = meshgrid(l,w);
h = 0.5*(W.^2).*(L.^2) + (c/2)*((L-x_l).^2+(W-x_w).^2);
figure;
contour(L,W,log10(h),15,'ShowText','off');
axis square;
end

function [f,g1,g2,H1,H2,H12] = f_eval_wl_matrix (L,W,X_l,X_w,rho,gamma)
	c = 1/(rho*gamma);
	f = 0.5*(L.^2).*(W.^2) + (c/2)*((L-X_l).^2+(W-X_w).^2);
	if nargout>=2
		g1 = zeros(size(L));
		g2 = zeros(size(L));
		g1 = (W.^2).*L + c*(L-X_l);
		g2 = (L.^2).*W + c*(W-X_w);
	end
	if nargout>=4
		H1 = W.^2+c;
		H2 = L.^2+c;
		H12 = 2*(L.*W);
	end
end


function [f,g,H] = f_eval_wl(y,x_l,x_w,rho,gamma)

l = y(1);
w = y(2);
c = 1/(rho*gamma);
f = 0.5*(l^2).*(w^2) + (c/2)*((l-x_l)^2+(w-x_w)^2);

g = zeros(2,1);
g(1,1) = w^2*l + c*(l-x_l);
g(2,1) = l^2*w + c*(w-x_w);

H = zeros(2,2);
H(1,1) = w^2+c;
H(2,2) = l^2+c;
H(1,2) = 2*l*w;
H(2,1) = 2*l*w;
end