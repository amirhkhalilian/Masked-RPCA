function [L,W,f,iters,grad_norm] = newtmeth_multisolver(fun,L0, W0, X_l, X_w, tol, maxit, c,modif_mode)

L = L0;
W = W0;
[f,g1,g2,H1,H2,H12] = fun(L,W);
grad_norm = sqrt((g1.^2)+(g2.^2));

for iters = 1:maxit
    work_ind = find(grad_norm>tol);
    work_size = length(work_ind);
    if isempty(work_ind)
        break;
    end
    % if iters<=50
    %     L(work_ind) = (c*X_l(work_ind))./(W(work_ind).^2+c);
    %     W(work_ind) = (c*X_w(work_ind))./(L(work_ind).^2+c);
    %     continue;
    % end
    % diag loading for problems with negative lambda
    lamda_min = 0.5*((L(work_ind)).^2+(W(work_ind)).^2+2*c)...
              -0.5*sqrt(((L(work_ind)).^2+(W(work_ind)).^2).^2+12*((W(work_ind)).^2).*((L(work_ind)).^2));    
    ind_lam_neg = find(lamda_min<0);
    switch modif_mode
        case 'AddAbs'
            H1(work_ind(ind_lam_neg)) = H1(work_ind(ind_lam_neg))+abs(lamda_min(ind_lam_neg));
            H2(work_ind(ind_lam_neg)) = H2(work_ind(ind_lam_neg))+abs(lamda_min(ind_lam_neg));
        case 'reverse'
            lamda_max = 0.5*((L(work_ind)).^2+(W(work_ind)).^2+2*c)...
              +0.5*sqrt(((L(work_ind)).^2+(W(work_ind)).^2).^2+12*((W(work_ind)).^2).*((L(work_ind)).^2));
            H1(work_ind(ind_lam_neg)) = lamda_max(ind_lam_neg);
            H2(work_ind(ind_lam_neg)) = lamda_max(ind_lam_neg);
            H12(work_ind(ind_lam_neg)) = 0;
        case 'significant'
            % calculate eigvals
            lamda_max = 0.5*((L(work_ind)).^2+(W(work_ind)).^2+2*c)...
              +0.5*sqrt(((L(work_ind)).^2+(W(work_ind)).^2).^2+12*((W(work_ind)).^2).*((L(work_ind)).^2));

            ind_d_h = find(lamda_min<0 & abs(H12(work_ind))<=1e-10);
            H1(work_ind(ind_d_h)) = lamda_max(ind_d_h);
            H2(work_ind(ind_d_h)) = abs(lamda_min(ind_d_h));
            H12(work_ind(ind_d_h)) = 0;
            
            ind_nd_h = find(lamda_min<0 & abs(H12(work_ind))>1e-10);
            V12 = lamda_max(ind_nd_h)-H1(work_ind(ind_nd_h));
            V22 = lamda_min(ind_nd_h)-H1(work_ind(ind_nd_h));
            V11 = H12(work_ind(ind_nd_h));
            V21 = H12(work_ind(ind_nd_h));
            n1 = sqrt(V11.^2+V21.^2);
            n2 = sqrt(V12.^2+V22.^2);
            V12 = V12./n2;
            V22 = V22./n2;
            V11 = V11./n1;
            V21 = V21./n1;
            lam_mod1 = lamda_max(ind_nd_h);
            lam_mod2 = abs(lamda_min(ind_nd_h));
            H1(work_ind(ind_nd_h)) = lam_mod1.*V11.^2+lam_mod2.*V12.^2;
            H2(work_ind(ind_nd_h)) = lam_mod1.*V21.^2+lam_mod2.*V22.^2;
            H12(work_ind(ind_nd_h)) = lam_mod1.*V11.*V21+lam_mod2.*V12.*V22;
        otherwise
            error('modification mode invalid!');
    end

    % Calculate the ascend direction as (H+tau*I)^{-1}*g with "\" operator
    % and the factorization result.
    det_hess = H1(work_ind).*H2(work_ind)-H12(work_ind).*H12(work_ind);
    dir1 = (1./det_hess).*(H2(work_ind).*g1(work_ind)-H12(work_ind).*g2(work_ind));
    dir2 = (1./det_hess).*(-H12(work_ind).*g1(work_ind)+H1(work_ind).*g2(work_ind));

    % Calculate t using BTLS. direction is -dir
    t = BTLS_multisolver(fun, L(work_ind), W(work_ind), X_l(work_ind), X_w(work_ind), c, -dir1, -dir2, 0.25, 0.5);

    % update 
    L(work_ind) = L(work_ind) - t.*dir1;
    W(work_ind) = W(work_ind) - t.*dir2;

    % calculate next values
    [f,g1,g2,H1,H2,H12] = fun(L,W);
    grad_norm = sqrt((g1.^2)+(g2.^2));

    % print 
    fprintf('i:%2d|mt:%2.4f|Mt:%2.4f|mg:%1.2e|Mg:%1.2e|w_size:%2d\n',iters,min(t),max(t),min(min(grad_norm)),max(max(grad_norm)),work_size);

end
end

function [t] = BTLS_multisolver(fun, L, W, X_l, X_w, c, dir1, dir2, alpha, rho)
    % check inputs
    if (alpha>1 || alpha<0)
        error('Invalid input c');
    end
    if (rho>1 || rho<0)
        error('Invalid input rho');
    end

    t = ones(size(dir1));

    while true
        % [f,g1,g2] = fun(L,W);
        % [f_next] = fun(L+t.*dir1,W+t.*dir2);
        f = 0.5*(L.^2).*(W.^2) + (c/2)*((L-X_l).^2+(W-X_w).^2);
        g1 = (W.^2).*L + c*(L-X_l);
        g2 = (L.^2).*W + c*(W-X_w);
        L_next = L+t.*dir1;
        W_next = W+t.*dir2;
        f_next = 0.5*(L_next.^2).*(W_next.^2) + (c/2)*((L_next-X_l).^2+(W_next-X_w).^2);
        ind_update_t = find(f_next>f+alpha*t.*(dir1.*g1+dir2.*g2));
        if isempty(ind_update_t)
            break;
        else
            t(ind_update_t) = rho*t(ind_update_t);
        end
    end

end



