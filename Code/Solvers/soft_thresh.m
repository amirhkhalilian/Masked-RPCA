function [Y_out] = soft_thresh(Y,tau)
	Temp = abs(Y)-tau;
	Temp = (Temp+abs(Temp))/2;
	Y_out = sign(Y).*Temp;
end