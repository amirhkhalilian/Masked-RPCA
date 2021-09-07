function [RE,SP,FPR,FNR,PWC,PR,F1] = CalculateStates(GT,Mask)
% This function calculates the statistics suggested by the CDnet paper. 

GT = GT(:);
Mask = Mask(:);

if length(GT)~=length(Mask)
	error('sizes do not mach');
end

TP = length(find(GT==1 & Mask==1));
TN = length(find(GT==0 & Mask==0));
FP = length(find(GT==0 & Mask==1));
FN = length(find(GT==1 & Mask==0));

RE = TP/(TP+FN);
SP = TN/(TN+FP);
FPR = FP/(FP+TN);
FNR = FN/(TN+FP);
PWC = 100*(FN+FP)/(TP+TN+FP+FN);
PR = TP/(TP+FP);
F1 = 2*(PR*RE)/(PR+RE);

end
