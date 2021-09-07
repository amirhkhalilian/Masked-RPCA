function  [X3D,W3D,X,W] = dataloader(dataset_name,mode,seque,varargin)
% function to loads the data from the CDnet 
if nargin<2
	fprintf('No data loading mode was provided to dataloader!\n');
	fprintf('useing auto generated directory!\n');
	mode = 'autodir';
	seque = [810,880,3];
end
switch lower(mode)
	case 'autodir'
		rootdir = pwd;
		datadir = [rootdir,filesep,'Data',filesep,dataset_name,filesep];
		inputdir = [datadir,'input',filesep];
		groundtruthdir = [datadir,'groundtruth',filesep];
	case 'datadir'
		inputdir = varargin{1};
		groundtruthdir = varargin{2};
	otherwise
		error('Invalid data loader mode!');
end
index=1;
for i = seque(1):seque(3):seque(2)
	in_name = sprintf('in%06d.jpg',i);
	gt_name = sprintf('gt%06d.png',i);
	X3D(:,:,index) = rgb2gray(imread([inputdir,in_name]));
	W3D(:,:,index) = imread([groundtruthdir,gt_name]);
	index = index+1;
end
X3D = double(X3D)/255; % bring the renge to [0,1]
W3D = double(W3D)/255; % bring the renge to [0,1]
if nargout>2
	[m,n,k] = size(X3D);
	X = reshape(X3D,m*n,k);
	W = reshape(W3D,m*n,k);
end

end