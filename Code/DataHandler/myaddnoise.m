function [Y3D] = myaddnoise(X3D,sigma)
[m,n,k] = size(X3D);
rng(0);
Y3D = X3D + sqrt(sigma)*randn(size(X3D));
end