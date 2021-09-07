function [Y_out] = svd_thresholding(Y,tau)
	% The slow function for SVT
	% This function computes the entire SVD and thresholds
	[U,S,V] = svd(Y,'econ');
	% S = sparse(diag(max(diag(S)-tau,0)));
	S = diag(S);
	S = sparse(diag(S(S>tau)-tau));
	Y_out = U(:,1:size(S,1)) * (S * V(:,1:size(S,1))'); 
end