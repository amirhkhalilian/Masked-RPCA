function [L,S,options] = RPCA_solver_my(X,options)

% Solver for ORIGINAL RPCA as following formulation:
% minimize |L|_{*} + lam |S|_1
% subject to:
% X = L + S

if nargin<2
	fprintf('Default initialization of parameters are set!\n');
	S = zeros(size(X));
	L = zeros(size(X));
% 	rho = .1;
	rho = numel(X)/(4*sum(abs(X(:))));
	maxiter = 100;
	tol = 1e-6;
	verbose = true;
else
	if isfield(options,'S')
		S = options.S;
	else
		S = zeros(size(X));
	end
	if isfield(options,'L')
		L = options.L;
	else
		L = zeros(size(X));
	end
	if isfield(options,'rho')
		rho = options.rho;
	else
		rho = numel(X)/(4*sum(abs(X(:))));
	end
	if isfield(options,'maxiter')
		maxiter = options.maxiter;
	else
		maxiter = 1000;
	end
	if isfield(options,'tol')
		tol = options.tol;
	else
		tol = 1e-6;
	end
	if isfield(options,'verbose')
		verbose = options.verbose;
	else
		verbose = true;
	end
end

% one line function for printing 
print_note = @(Str) fprintf([repmat('-',1,80),'\n',Str,'\n',repmat('-',1,80),'\n']);
% initialization of the Dual Variables
U = zeros(size(X));
lam = .08/sqrt(max(size(X)));

% set flag to report convergence
flag_tol_reached = false;
print_note(sprintf('Starting algorithm with %4d maximum iterations',maxiter));

% start timing
tic
% print initial values of the objective function
if verbose
	norm_nuc_L(1) = sum(svd(L,'econ'));
	norm_1_S(1) = sum(sum(abs(S)));
	gap(1) = norm(X - (L+S),'fro');
	obj_total(1) = norm_nuc_L(1) + lam*norm_1_S(1);
	fprintf([repmat('-',1,80),'\n']);
	fprintf('i:%3d|obj:%1.2e|gap:%1.2e|nucL:%1.2e|nS:%1.2e|\n',...
				0,obj_total(1),gap(1),norm_nuc_L(1),norm_1_S(1));
end

for iter = 1:maxiter
	% Update L
	L = svd_thresholding(X-S+U/rho,1/rho);
	% Update S
	S = soft_thresh(X-L+U/rho,lam/rho);
	% Update U
	U = U + rho*(X-L-S);
	% convergence
	gap(iter+1) = norm(X - (L+S),'fro');
	if verbose
		norm_nuc_L(iter+1) = sum(svd(L,'econ'));
		norm_1_S(iter+1) = sum(sum(abs(S)));
		obj_total(iter+1) = norm_nuc_L(iter+1) + lam*norm_1_S(iter+1);
		fprintf('i:%3d|obj:%1.2e|gap:%1.2e|nucL:%1.2e|nS:%1.2e|\n',...
					iter,obj_total(iter+1),gap(iter+1),norm_nuc_L(iter+1),norm_1_S(iter+1));
	end
	if gap(iter+1)<=tol
		break;
	end
end

% record total time
time_total = toc;

% return the convergance properties in the options struct
options.time_total = time_total;
options.norm_nuc_L = norm_nuc_L;
options.norm_1_S = norm_1_S;
options.gap = gap;
options.obj_total = obj_total;
options.flag_tol_reached = flag_tol_reached;
options.iterations = iter;

% print a summary of the convergance results
if verbose
	print_note('Summary of the results:');
	if flag_tol_reached
		print_note(sprintf('tolerance reached with %4d iterations',iter));
	else
		print_note(sprintf('terminated with maximum iteration without reaching the tolerance'));
	end
	fprintf('total time elapsed: %3.3f\n',time_total);
	fprintf('final objective value: %1.5e\n',obj_total(end));
	fprintf('Duality gap || X - (L+S) ||: %1.5e\n',gap(end));
	fprintf('|| L ||* : %1.5e\n',norm_nuc_L(end));
	fprintf('|| S ||1 : %1.5e\n',norm_1_S(end));
	print_note('End of Summary!');
end
end