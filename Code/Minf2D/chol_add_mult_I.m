function [R,tau] = chol_add_mult_I(H,beta)
% this function implements the Cholesky with added multiple of the identity.
% as in N&W algorithm 3.3.

% Inputs:
% H: Matrix to be factorized (Hessian)
% beta: 

% output

[m,n] = size(H);
I = speye(m,n);

if min(diag(H))>0
	tau = 0;
else
	tau = -min(diag(H))+beta;
end

while true;
	[R,p] = chol(H+tau*I);
	if ~p
		break;
	else
		tau = max([2*tau, beta]);
	end
end

end