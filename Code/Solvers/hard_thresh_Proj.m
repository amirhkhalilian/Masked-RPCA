function [Y_out] = hard_thresh_Proj(X,lambda)
	% Y_out = X.*(abs(X)>tau);
	Y_out = zeros(size(X));
	Ind_1 = find(X<=1 & X>=0 & X>=sqrt(2*lambda));
	Y_out(Ind_1) = X(Ind_1);
	Ind_2 = find(X>=1 & X>=(lambda+0.5));
	Y_out(Ind_1) = 1;
end