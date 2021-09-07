function [L,W] = stupid_multisolver(fun, L0, W0, X_l, X_w, tol, maxit, c)
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