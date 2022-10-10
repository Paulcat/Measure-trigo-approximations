function P = proxy_nonlinear(c,n,options)
%PROXY_NONLINEAR compute a  polynomial proxy from moment matrix
%   P = PROXY_NONLINEAR(M,N,OPTIONS) compute a polynomial from the singular
%   vectors of the moment matrix M.
%
%   N is the order of the moment matrix.
%
%   OPTIONS:
%		- sigma (default: 1e-5):
%			- if integer: exact rank of moment matrix (or at least upper bound)
%			- if not: tolerance for numerical rank
%		- grid_size (default: 101)


d = numel(n);
N = prod(n+1);
M = prod(options.grid_size);

c = reshape(c,[2*n+1]);
c = ifftshift(c); % 0 frequency in the middle


% load options
sigma = getoptions(options,'sigma',1e-5);
mgrid = getoptions(options,'grid_size',101*ones(1,d));
draw  = getoptions(options,'draw',0); 


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

% compute proxy (naive implementation)
%Xf = compute_grid(n,'spectral');
%G  = compute_grid(mgrid,'spatial'); L = prod(mgrid);
%F = exp(-2i*pi*sum( reshape(Xf,[N,1,d]).*reshape(G,[1,L,d]),3));
%P = sum(abs(F'*U).^2,2);
%P = reshape(P,[mgrid,1]);

% (fft implementation)
U  = reshape(U,[n+1,r]);
Uf = M * ifft2( padarray(U,mgrid-n-1,'post') ); % padding along all dimension but last
%TODO: implementation for arbitrary dimension. This only works for d=2
P  = sum(abs(Uf).^2,d+1);
P  = P/sum(P(:));


if draw
	G = compute_grid(mgrid,'spatial');
	clf;
	if d==1
		plot(G,P,'linewidth',3);
	elseif d==2
		surf(G(:,:,1),G(:,:,2),P,'linestyle','none');
		view(2);
		drawnow;
	else
		warning('Cannot display for dimensions higher thna 2');
	end
end
