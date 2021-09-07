clc; clear; close all;

seque = [1250,1630,1];
gtseque = [1499,1515,1523,1547,1548,1553,1554,1559,1575,1577,1594,1597,1601,1605,1615,1616,1619,1620,1621,1624];
rootdir = pwd;
datadir = [rootdir,filesep,'Data',filesep,'watersurface',filesep];
inputdir = [datadir];
index=1;
for i = seque(1):seque(3):seque(2)
	in_name = sprintf('WaterSurface%4d.bmp',i);
	% gt_name = sprintf('gt%06d.png',i);
	X3D(:,:,index) = rgb2gray(imread([inputdir,in_name]));
	% W3D(:,:,index) = imread([groundtruthdir,gt_name]);
	index = index+1;
end
for i = 1:length(gtseque)
	gt_name = sprintf('gt_new_WaterSurface%4d.bmp',gtseque(i));
	W3D(:,:,i) = rgb2gray(imread([inputdir,gt_name]));
end	
X3D = double(X3D)/255; % bring the renge to [0,1]
W3D = double(W3D)/255; % bring the renge to [0,1]

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

options.lam = 20/(max(size(X))); % good parameters
options.rho =  1e-1;
options.gam = 30;
options.maxit =200;
options.adaptive = true;

[Yl,Yw,Zl,Zw,Ul,Uw,options] = Lnuc_L0_DR_solver(X,L0,W0,options);

figure;
subplot(131); plot(options.D); title('DR merit function');
subplot(132); plot(options.objY); hold on; plot(options.objZ); title('objective value');
subplot(133); plot(options.diff_ZY); title('Z-Y');


Visual_view3d(reshape(Zw,m,n,k),reshape(Yw,m,n,k));
Visual_view3d(reshape(Zl,m,n,k),reshape(Yl,m,n,k));


% save('watersurface_Res.mat');