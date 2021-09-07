function plot3D_showpair(names,varargin)

number_in = length(varargin);
figure;
ha = tight_subplot(1,number_in,0,[0,0.03],0);
for i = 1:number_in
	[m,n,k] = size(varargin{i});
	X2Dcolumn = [];
    X_GT = [];
	for frame = 1:k
		X2Dcolumn = cat(1,X2Dcolumn,squeeze(varargin{i}(:,:,frame)));
        X_GT = cat(1,X_GT,squeeze(varargin{1}(:,:,frame)));
	end
	axes(ha(i));
	imshowpair(X_GT,X2Dcolumn); colormap(gray);
	axis off;
    title(names{i});
end
end