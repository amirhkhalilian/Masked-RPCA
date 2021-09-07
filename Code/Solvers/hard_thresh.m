function [Y_out] = hard_thresh(X,tau)
	Y_out = X.*(abs(X)>tau);
end