clc; clear; close all;

mode = 'autodir';
seque = [630,740,1];
[X3D,W3D] = dataloader('pedestrians',mode,seque);
[m,n,k] = size(X3D);
X = reshape(X3D,m*n,k);
W_true = reshape(W3D,m*n,k);


[m,n,k] = size(X3D);
X = reshape(X3D,m*n,k);
L0 = median(X3D,3);
S0 = X3D - repmat(L0,[1,1,k]);
L0 = repmat(L0(:),[1,k]);
W0 = reshape(abs(S0),m*n,k);

% options.lam = 1/(max(size(X))); % good parameters
% options.rho = 1e-2;
% options.gam = 10;
% options.maxit =100;

options.lam = 1/(max(size(X))); % good parameters
options.rho =  1e-1;
options.gam = 50;
options.maxit =20;
options.adaptive = true;

[Yl,Yw,Zl,Zw,Ul,Uw,options] = Lnuc_L0_DR_solver(X,L0,W0,options);

figure;
subplot(131); plot(options.D); title('DR merit function');
subplot(132); plot(options.objY); hold on; plot(options.objZ); title('objective value');
subplot(133); plot(options.diff_ZY); title('Z-Y');


Visual_view3d(reshape(Zw,m,n,k),reshape(Yw,m,n,k));
Visual_view3d(reshape(Zl,m,n,k),reshape(Yl,m,n,k));