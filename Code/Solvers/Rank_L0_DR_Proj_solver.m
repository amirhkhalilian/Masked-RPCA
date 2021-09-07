function [Yl,Yw,Zl,Zw,Ul,Uw,options] = Rank_L0_DR_Proj_solver(X,L0,W0,options)
	% this function solves the minimization problem
	% ||L||_* + lam ||W||_0 + (rho/2) ||(W-1)o(L-X)||^2
	% using DR-method

	lam = options.lam;
	rho = options.rho;
	gam = options.gam;
	maxit = options.maxit;
	change_gam = options.adaptive;
	Ul = L0;
	Uw = W0;
	Yl = L0;
	Yw = W0;

	for iter = 1:maxit
		% minimize f(y) + 1/(2gam)||y-u||^2
		fun_1 = @(L,W)f_eval_wl_matrix(L,W,Ul,Uw,rho,gam);
		[Yl,Yw] = multisolver_GS(fun_1, Yl-X, Yw-1, Ul-X, Uw-1, 1e-5, 100, 1/(gam*rho));
		Yl = Yl+X;
		Yw = Yw+1;
		Yw(Yw<0) = 0;
		Yw(Yw>1) = 1;
		% minimize g(z) + 1/(2gam)||z-2y+u||^2
		Zl = svd_hardthresholding(2*Yl-Ul,gam);
		Zw = hard_thresh_Proj(2*Yw-Uw,gam*lam);
		% Zw = hard_thresh(2*Yw-Uw,sqrt(2*gam*lam));
		% Zw(Zw<0) = 0;
		% Zw(Zw>1) = 1;

		% calculate Douglas Rachford merit function
		% nnuc_Zl = sum(svd(Zl,'econ'));
		% n0_Zw = nnz(Zw);
		% g(iter) = nnuc_Zl + lam* n0_Zw;
		% f(iter) = (rho/2)*norm((1-Yw).*(Yl-X),'fro');
		% D(iter) = f(iter)+g(iter) -(1/(2*gam))*(norm(Yl-Zl,'fro')^2+norm(Yw-Zw,'fro')^2) + ...
		% 		  (1/gam)*(sum((Ul(:)-Yl(:)).*(Zl(:)-Yl(:)))+sum((Uw(:)-Yw(:)).*(Zw(:)-Yw(:))));

		nnuc_Zl = sum(svd(Zl,'econ'));
		n0_Zw = nnz(Zw);
		nnuc_Yl = sum(svd(Yl,'econ'));
		n0_Yw = nnz(Yw);
		gZ(iter) = nnuc_Zl + lam* n0_Zw;
		fY(iter) = (rho/2)*norm((1-Yw).*(Yl-X),'fro');
		gY(iter) = nnuc_Yl + lam* n0_Yw;
		fZ(iter) = (rho/2)*norm((1-Zw).*(Zl-X),'fro');
		D(iter) = fY(iter)+gZ(iter) -(1/(2*gam))*(norm(Yl-Zl,'fro')^2+norm(Yw-Zw,'fro')^2) + ...
				  (1/gam)*(sum((Ul(:)-Yl(:)).*(Zl(:)-Yl(:)))+sum((Uw(:)-Yw(:)).*(Zw(:)-Yw(:))));

		objY(iter) = gY(iter)+(rho/2)*fY(iter);
		objZ(iter) = gZ(iter)+(rho/2)*fZ(iter);
		diff_ZY(iter) = sqrt(norm(Yl(:)-Zl(:))^2+norm(Yw(:)-Zw(:))^2);
		fprintf('i:%2d|D:%1.4e|f:%1.2e|g:%1.2e|\n',iter,D(iter),fY(iter),gZ(iter));

		% update u		
		Ul = Ul+(Zl-Yl);
		Uw = Uw+(Zw-Yw);
		if change_gam
			gam = max([10,gam*0.99]);
		else
			gam = gam;
		end
		
	end
	options.D = D;
	options.fY = fY;
	options.gY = gY;
	options.objY = objY;
	options.fZ = fZ;
	options.gZ = gZ;
	options.objZ = objZ;
	options.diff_ZY = diff_ZY;

end

function [L,W] = multisolver_GS(fun, L0, W0, X_l, X_w, tol, maxit, c)
check_hess = true;
L = L0;
W = W0;
[f,g1,g2,H1,H2,H12] = fun(L,W);
det_hess = H1.*H2-H12.^2;
grad_norm = sqrt((g1.^2)+(g2.^2));

for iters = 1:maxit
    work_ind = find(grad_norm>tol);
    work_size = length(work_ind);
    if isempty(work_ind)
        break;
    end
  L(work_ind) = (c*X_l(work_ind))./(W(work_ind).^2+c);
  W(work_ind) = (c*X_w(work_ind))./(L(work_ind).^2+c);
  g1(work_ind) = (W(work_ind).^2).*L(work_ind) + c*(L(work_ind)-X_l(work_ind));
  g2(work_ind) = (L(work_ind).^2).*W(work_ind) + c*(W(work_ind)-X_w(work_ind));
  grad_norm(work_ind) = sqrt((g1(work_ind).^2)+(g2(work_ind).^2));
  if check_hess
  	H1(work_ind) = W(work_ind).^2+c;
	H2(work_ind) = L(work_ind).^2+c;
	H12(work_ind) = 2*(L(work_ind).*W(work_ind));
	det_hess(work_ind) = H1(work_ind).*H2(work_ind)-H12(work_ind).^2; 
  end
  fprintf('i:%3d|mg:%1.2e|Mg:%1.2e|mH:%+1.2e|MH:%+1.2e|nH:%4d|w_size:%2d\n',...
  		iters,min(min(grad_norm)),max(max(grad_norm)),min(min(det_hess)),max(max(det_hess)),length(find(det_hess<0)),work_size);
end
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
