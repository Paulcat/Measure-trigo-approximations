function P = proxy_nonlinear(c,n,options)
%PROXY_NONLINEAR compute a  polynomial proxy from moment matrix
%   P = PROXY_NONLINEAR(M,N,OPTIONS) compute a polynomial from the singular
%   vectors of the moment matrix M.
%
%   N is the order of the moment matrix.
%
%   Supported OPTIONS:
%	- order (monomial ordering)  : colex(*) | lex
%	- svd_tol (tolerance for svd): 1e-8(*)
%	- grid_size


d = numel(n);
N = prod(n+1);
M = prod(options.grid_size);


% load options
order = getoptions(options,'order','colex');
sigma = getoptions(options,'rank',1e-5);
mgrid = getoptions(options,'grid_size',101*ones(1,d));


% make sure grid is sufficiently fine wrt cutoff frequency
gc = ceil((mgrid-1)/2);
if sum(gc<n+1)
	cut  = sprintf('%d', n+1);
	grid = sprintf('%d', gc);
	error('Grid coarser than cutoff frequency: grid=%s, fc=%s', grid, cut);
end

% compute SVD
%[U,S] = svd(M); S = diag(S);
if mod(sigma,1)==0 % user-specified rank
	r = min(sigma,N); %TODO: inexact value of rank for curve?
	[U,S] = svds(@(v,tflag)Tprod(c,reshape(v,[n+1,1]),tflag),[N,N],r);
	S = diag(S);
	r = find(S/max(S) < 1e-8,1) -1; % further refinement
	if isempty(r)
		r = length(S);
	end
	U = U(:,1:r);
else % numerical rank
	[U,S] = svds(@(v,tflag)Tprod(c,reshape(v,[n+1,1]),tflag),[N,N],N);
	S = diag(S);
	r     = find(S/max(S) < sigma,1) - 1; % numerical rank
	if isempty(r) %TODO:HACK
		r  = length(S);
	end
	U    = U(:,1:r);
end

% compute proxy
%Xf = compute_grid(n,'spectral');
%G  = compute_grid(mgrid,'spatial'); L = prod(mgrid);
%F = exp(-2i*pi*sum( reshape(Xf,[N,1,d]).*reshape(G,[1,L,d]),3));
%P = sum(abs(F'*U).^2,2);
%P = reshape(P,[mgrid,1]);

% with fft
U  = reshape(U,[n+1,r]);
Uf = M * ifft2( padarray(U,mgrid-n-1,'post') ); % padding along all dimension but last
%TODO: implementation for arbitrary dimension. This only works for d=2
P  = sum(abs(Uf).^2,d+1);
P  = P/sum(P(:));

end
