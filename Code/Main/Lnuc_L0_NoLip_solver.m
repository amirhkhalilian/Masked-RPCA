function [L,W,options] = Lnuc_L0_NoLip_solver(X,L0,W0,options)
	% this function solves the minimization problem
	% ||L||_* + lam ||W||_0 + (rho/2) ||(W-1)o(L-X)||^2
	% using the NoLip method.
	lam = options.lam;
	rho = options.rho;
	gam = options.gam;
	maxit = options.maxit;
	L = L0;
	W = W0;
	L_prev = L0;
	W_prev = W0;

	gam_min = sqrt((rho^2)/4+3*rho) + rho;
	if gam<=gam_min
		fprintf('the gamma is less than gam_min = %1.2f\n',gam_min);
		% return;
	end

	for iter = 1:maxit
		Psi_L = ((W-1).^2).*(L-X);
		Ak = X - gam*Psi_L + (L-X);
		Psi_W = (W-1).*((L-X).^2);
		Bk = 1 - gam*Psi_W + (W-1);
		L = svd_thresholding(Ak,gam);
		W = hard_thresh(Bk,sqrt(2*gam*lam));
		W(W<0) = 0;
		W(W>1) = 1;
		L_nuc = sum(svd(L,'econ'));
		n0_W = nnz(W);
		obj(iter) = L_nuc + lam*n0_W;
		L_change(iter) = norm(L-L_prev,'fro');
		W_change(iter) = norm(W-W_prev,'fro');
		fprintf('i:%d|obj:%1.4e|cL:%1.2e|cW:%1.2e\n',iter, obj(iter),L_change(iter),W_change(iter));
		L_prev = L;
		W_prev = W;
	end

	options.obj = obj;
	options.L_change = L_change;
	options.W_change = W_change;

end