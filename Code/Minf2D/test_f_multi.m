clc; clear; close all;

x_l_1 = 0.1;
x_w_1 = 0.1;
gam = 2;
rho = 2;
maxL = 4;
maxW = 4;

x_l_2 = 1;
x_w_2 = 1;
gam = 2;
rho = 2;
maxL = 4;
maxW = 4;

x_l_3 = 10;
x_w_3 = 10;
gam = 2;
rho = 2;
maxL = 15;
maxW = 15;

x_l_4 = 3;
x_w_4 = 10;
gam = 2;
rho = 2;
maxL = 15;
maxW = 15;

X_l = [x_l_1;x_l_2;x_l_3;x_l_4];
X_w = [x_w_1;x_w_2;x_w_3;x_w_4];

X_l = randn(1000,1000);
X_w = randn(1000,1000);


minC = -0.5*(maxL^2+maxW^2)+sqrt(3*(maxL^2)*(maxW^2)+0.25*(maxL^2+maxW^2)^2);

fprintf('minimmum C value is: %+1.3f\n chosen c is: %+1.3f\n',minC,1/(rho*gam));

fun = @(L,W)f_eval_wl_matrix(L,W,X_l,X_w,rho,gam);

% [z_nmeth_1,f_nmeth,iters,norm_grad] = newtmeth(fun, [2;-4], 1e-9, 100);
% [z_nmeth_2,f_nmeth,iters,norm_grad] = newtmeth_eigmod(fun, [2;-4], 1e-9, 100,1/(gam*rho));
% [z_nmeth_3,f_nmeth,iters,norm_grad] = newtmeth_eigabs(fun, [2;-4], 1e-9, 100);
% [f,g,H] = f_eval_wl(z_nmeth_1(:,end),x_l,x_w,rho,gam);
% eig(H)

% [L,W,f,iters,grad_norm] = newtmeth_multisolver_fast(fun, [2;2;2;2], [-4;-4;-4;-4], X_l,X_w, 1e-9, 100, 1/(gam*rho),'significant');
[L,W,f,iters,grad_norm] = newtmeth_multisolver_fast(fun, randn(1000,1000),randn(1000,1000), X_l,X_w, 1e-9, 1000, 1/(gam*rho),'reverse');
% [L,W] = stupid_multisolver(fun, randn(1000,1000),randn(1000,1000), X_l,X_w, 1e-9, 1000, 1/(gam*rho));

% plot_2D_contour(maxL,maxW,x_l,x_w,rho,gam)
% hold on;
% plot(z_nmeth_1(1,:),z_nmeth_1(2,:),'r-*','MarkerSize',10);
% plot(z_nmeth_2(1,:),z_nmeth_2(2,:),'b-*','MarkerSize',10);
% plot(z_nmeth_3(1,:),z_nmeth_3(2,:),'m-*','MarkerSize',10);

% plot(z_nmeth(1,:),z_nmeth(2,:),'-*')

function plot_2D_contour(maxL,maxW,x_l,x_w,rho,gamma)
c = 1/(rho*gamma);
l = -maxL:0.01:maxL;
w = -maxW:0.01:maxW;
[L,W] = meshgrid(l,w);
h = 0.5*(W.^2).*(L.^2) + (c/2)*((L-x_l).^2+(W-x_w).^2);
figure;
contour(L,W,log10(h),15,'ShowText','on');
axis square;
end

function [f,g1,g2,H1,H2,H12] = f_eval_wl_matrix(L,W,X_l,X_w,rho,gamma)
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
