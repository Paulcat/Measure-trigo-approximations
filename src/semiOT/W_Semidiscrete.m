function [OT,w] = W_Semidiscrete(G,f,x,a,options)
%W_SEMIDISCRETE Compute semi-discrete optimal transport
%   OT = WASSERSTEIN_DIST(G,F,X,A,OPT) computes the Wasserstein distance
%   between the density F, evaluated on the grid G, and the discrete
%   measure supported on X, with amplitudes A.
%
%   Main options:
%		  maxit (100) 			-- maximum number of iterations for descent
%       method ('l-bfgs') 	-- solver
%		  multiscale (1) 		-- trigger multiscale implementation
%		  norm (1) 				-- choice of norm for OT
%		  parallelize (0) 	-- trigger parallelized implementation
%		  wasserstein (1) 	-- choice of wasserstein distance

%disp('Computing semidiscrete OT...');

d = size(x,2);
s = size(x,1);
m = size(G,1);


%TODO: parallelization


% safety nets
if d ~= size(G,2)
    error(['Dimension mismatch between density (', int2str(size(G,2)), ') '...
        'and discrete (', int2str(d), ')']);
end

if m ~= numel(f)
	error(['Mismatch in the number of points for density: amplitude vector is ', ...
		int2str(numel(f)), ' but the grid has ' int2str(size(G,1)), ' points']);
end


% load options
[verbose,verboseI,method,maxit,tau,progtol,p,q,nbloc,L] = ...
	processOptions(options);


% helpers
flat = @(x) x(:);
vout = @(varargin) deal(varargin{1:nargout});


% multiscale decomposition of (x,a)
S		  = ceil(s./(2.^(0:L)));
tic;
[X,A,P] = multiscale_decomposition(x,a,L,S,q,p);
time = toc;

if verboseI
	fprintf('Quantization in %i levels took %.2d seconds\n',L,time);
	fprintf('Number of points are: [');
	fprintf('%g, ', S(1:end-1)); fprintf('%i]\n',S(end));
end

% dual potential
w = zeros(S(end),1);

l = L+1;
while l >= 1
%for l=L-1:-1:0 % multiscale levels, starting from avant-dernier
	% initialize
	s_l = S(l);
	x_l = X{l};
	a_l = A{l};

	if ~nbloc
		% compute cost matrix
		D = dist_matrix(G,x_l,p).^q;

		% semidiscrete transport functionals
		Ctrf = @(w) min(D-w(:)',[],2); % (semidiscrete) c-transform
		obj  = @(w) w(:)'*a_l(:) + Ctrf(w)'*f(:); % objective
		%
		Ids  = reshape(1:s_l,[1 1 s_l]);
		r 	  = @(w) flat(sum( (Laguerre_map(D,w)==Ids) .* f(:) )); % measure of Laguerre cells
		grad = @(w) a_l(:) - r(w); % gradient
		%
		fg  = @(w) vout(-obj(w),-grad(w));
		%fg = @(w) funSD(G,f,x_l,a_l,p,q,w,3);

	else
		f = f(:);
		% solution 1: unequal blocs, use cell arrays
		B 		  = bloc_divide(m,nbloc);
		Gsliced = cell(nbloc,1);
		fsliced = cell(nbloc,1);
		Dsliced = cell(nbloc,1);
		for i=1:nbloc
			bi 		  = B(i):(B(i+1)-1);
			Gsliced{i} = G(bi,:);
			fsliced{i} = f(bi);
			Gi = G(bi,:);
			Dsliced{i} = dist_matrix(Gi,x_l,p).^q;
		end

		% solution 2: enforce divisibility
		%if mod(m,nbloc)
		%	error(['please make sure that m (',int2str(m),') can be divided into n (' ...
		%		,int2str(nbloc),') equal parts']);
		%end
		%Gsliced = reshape(G',[d,m/nbloc,nbloc]);
		%Gsliced = permute(Gsliced,[2 1 3]);
		%fsliced = reshape(f,[m/nbloc,nbloc]);

		fg  = @(w) funSD_parallel(Gsliced,fsliced,x,a_l,p,q,w,nbloc);

	end

	switch method
		case 'g-ascent'
			OT = zeros(maxit,1);
        	for i=1:maxit
				w = w + tau*grad(w);
            OT(i) = obj(w).^(1/q);

				%m = sqrt(size(G,1));
				%t = (0:m-1)'/m;
				%I = Laguerre_map(D,w);
				%I = reshape(I,[m,m]);
				%imagesc(t,t,I);
				%contour(t,t,I,.5:s+.5,'r','Linewidth',3);
				%scatter(x_l(:,2),x_l(:,1),50,'r','filled');
				%xlim([0,1]), ylim([0,1]);
				%drawnow;
				%pause(1);
        	end

		case 'l-bfgs'
        	opt.MaxIter = maxit;
        	opt.display = 'full';
			opt.progTol = progtol; % tolerance for lack of progress
         [w,~,~,output] = lbfgs(fg,w,opt);
        	OT = (-output.trace.fval).^(1/q);
			mean_time = mean(output.trace.time);

			if verboseI
				fprintf('Total time: %.2fs\n',sum(output.trace.time));
			end

		otherwise
			error('Unknown method');
    end

	 % use new dual potential as next initialization
	 if ~isempty(P) && l>1 % not if this the last level
		 w = w(P{l-1});
	 end

	 l = l-1; % go to next scale

	 % display Laguerre cells
	 %if d==2 && display_Lag
	 if 0
		 m = sqrt(size(G,1));
		 t = (0:m-1)'/m;

		 % display Laguerre cells
		 clf, hold on;
		 I = Laguerre_map(D,w);
		 I = reshape(I,[m,m]);
		 imagesc(t,t,I);
		 contour(t,t,I,.5:s+.5,'r','Linewidth',3);
		 scatter(x_l(:,2),x_l(:,1),50,'r','filled');
		 %voronoi(x(:,2),x(:,1));
		 xlim([0,1]); ylim([0,1]);
		 drawnow;
		 pause(1);
	end
end

%disp('...end');
end


function [verbose,verboseI,method,maxit,tau,progtol,p,q,nbloc,L] = processOptions(options)

method  = getoptions(options,'method','l-bfgs'); % choice of descent algorithm
maxit   = getoptions(options,'maxit',200); % max iterations for descent
tau     = getoptions(options,'tau',.1); % step size for gradient descent
progtol = getoptions(options,'progtol',1e-9); % tolerance for descent wrt param changes
%
q       = getoptions(options,'wasserstein',1); % choice of Wasserstein distance
p       = getoptions(options,'norm',1); % choice of norm
%
nbloc   = getoptions(options,'parallelize',0); % number of blocs for parallel computing
L  	  = getoptions(options,'multiscale',0); % number of scales
%
%display_Lag = getoptions(options,'display_Lag',0); % display Laguerre cells
display = getoptions(options,'display','off'); % outstream: 'off' | 'reduced' | 'debug'

switch display
	case 'off'
		verbose  = 0;
		verboseI = 0;
		
	case 'reduced'
		verbose  = 1;
		verboseI = 0;

	case 'debug'
		verbose  = 1;
		verboseI = 1;
	
end

end


% obsolete code, for symmetrizing domain...
% Laguerre cells
% % symmetrize domain
% ts = linspace(-1,2,3*m);
% [Ys,Xs] = meshgrid(ts);
% Gs = [Xs(:),Ys(:)];
% %
% [s2,s1] = meshgrid(-1:1);
% shift = [s1(:),s2(:)];
% xs = reshape(x,[s,1,2]) + reshape(shift,[1,9,2]);
% xs = reshape(xs,[],2);
% as = repmat(a,9,1);
% %
% Ps = [P P P; P P P; P P P];

% display
% clf, hold on;
% scatter(xs(:,1),xs(:,2),40,'filled');
% plot([0 0],[-1 2],'k');
% plot([1 1],[-1 2],'k');
% plot([-1 2],[0 0],'k');
% plot([-1 2],[1 1],'k');
