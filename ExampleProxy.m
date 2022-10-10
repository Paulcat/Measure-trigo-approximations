% Code to reproduce the figures in [Catala, Hockmann, Kunis, PAMM 2022]

clear all
addpath('examples/','examples/data');
addpath('src/lbfgs', 'src/proxys', 'src/semiOT', 'src/toolbox');

if ~exist('results','dir')
	mkdir('results');
end


d = 2; % dimension



% ***** USER-DEFINED ***** %
% Measure to recover
type = 'd'; % d: discrete | z: curve | c: circle

% choice of proxy
proxy = 1; % 0: convolution | 1: signal polynomial
if ~proxy
	weights = 'F'; % F: Fejer | K2: Jackson | K3: Fejer^3
else
	weights = '-'; % (useful only for naming file when saving figures)
end
% ************************* %




% define the stepsize in n for the plots
if strcmp(weights, 'K3')
	step_n = 3; % take n of the form n = 3m only
	mycolor = 'blue'; % for plotting
elseif strcmp(weights, 'K2')
	step_n = 2; % only even n
	mycolor = 'blue';
elseif strcmp(weights, 'F')
	step_n = 2;
	mycolor = 'red';
elseif strcmp(weights, '-')
	step_n = 2;
	mycolor = 'red';
end




% ***** USER_DEFINED ***** %
% range of frequency to test
nmin = 1;
L	  = step_n * 15; % increasing frequency step_n by step_n
nmax = L + nmin -1;
% ************************ %




% instantiate measure
if strcmp(type,'d')
	s = 15; % sparsity
	load('discrete-measure.mat'); % use the same discrete measure every time
	mvec = @(kv) exp(-2i*pi*sum(reshape(kv,[],1,d) .* reshape(x,[1,s,d]),3))*a;

	% display
	clf, scatter(x(:,1),x(:,2),100*a,'filled');
	xlim([0,1]),ylim([0,1]);
else
	s 	  			= 1000; % sparsity. To use pre-computed data, choose 1000, 2000 or 3000
	samp 			= 'regular'; % sampling of the curve: 
	[x,a,mvec] 	= instantiate_measure(type,s,d,samp);

	% display
	clf, scatter(x(:,1),x(:,2),50,'.');
	xlim([0,1]),ylim([0,1]);
end


% discretization grid
mgrid = [100 100]; % grid size
Y 		= compute_grid(mgrid,'spatial'); % discretization grid
G		= reshape(Y,[prod(mgrid),d]);




% ***** OPTIONS FOR SEMIDISCRETE OPTIM ***** %
options_W.maxit = 100;
%options_W.parallelize = 0; %800
%options_W.multiscale = 0;
options_W.wasserstein = 1;
options_W.norm = options_W.wasserstein; % HACK
options_W.progtol = 1e-9;
% ****************************************** %


W 		= zeros(1,(nmax-nmin+1)/step_n);
% parpool(10);
% parfor n=1:(nmax-nmin+1) % it is possible to parallelize this loop
for n=1:(nmax-nmin+1)/step_n
	nvec = (step_n*n + (nmin-1))*ones(1,d); % e.g.  n = 2,4,6,...,250 
	N	  = prod(2*nvec+1);

	% set proxy options
	options = struct;
	options.grid_size = mgrid;
	options.kernel 	= weights;
	options.sigma 		= 2*sum(nvec); % (upper bound) only useful for computing svd
	options.draw 		= 1;

	% moment vector
	Xf = compute_grid(nvec,'spectral-sym');
	Xf = reshape(Xf,[N,d]);
	c  = mvec(Xf); % moment vector of the measure to recover

	if ~proxy % convolution polynomial
		P  = proxy_linear(c,nvec,options);
		if strcmp(type,'z') %TODO:HACK!
			P = fftshift(P);
		end
		OT = W_Semidiscrete(G,P,x,a,options_W);

	else % signal polynomial
		P  = proxy_nonlinear(c,nvec,options);
		if strcmp(type,'z') %TODO:HACK!
			P = fftshift(P);
		end
		% "unweighted" amplitudes
		a1 = ones(size(x,1),1);
		a1 = a1/norm(a1,1);
		OT = W_Semidiscrete(G,P,x,a1,options_W);
	end
	%OT = W_Semidiscrete(G,P,x,a,options_W);

	W(n) = OT(end);
end
fprintf('\n');
%delete(gcp('nocreate'));

filename = ['results/W', int2str(options_W.wasserstein), '_prox',...
	int2str(proxy),'-', weights,'_', type, '-', int2str(s), '_n', ...
	int2str(nmin), '-', int2str(nmax)];
fullname = [filename, '.csv'];
if isfile(fullname)
	csvwrite([filename,'-bis.csv'],W);
else
	csvwrite(fullname,W);
end


% display
nlist = step_n*nmin : step_n : nmax;
p1 = loglog(nlist,W,'linewidth',3,'color',mycolor);
hold on;
p2 = loglog(nlist,1./nlist,'--b','linewidth',2);
p3 = loglog(nlist,log(nlist)./nlist,'--r','linewidth',2);
p4 = loglog(nlist,1./sqrt(nlist),'-.r','linewidth',2);
p5 = loglog(nlist,1./(nlist).^(1/3),':r','linewidth',2);
legend([p1,p2,p3,p4,p5],...
	['$W_',int2str(options_W.wasserstein),'(',weights,', \mu_{',type '})$'],...
	'$1/n$','$\log(n)/n$','$1/\sqrt{n}$','$1/\sqrt[3]{n}$','interpreter','latex');
set(gca,'fontsize',25);