function P = proxy_linear(c,n,options)
%PROXY_LINEAR compute a polynomial proxy from moments
%   P = PROXY_LINEAR(C,N,OPTIONS) compute the convolution between a
%   measure and a band-limited kernel specified by OPTIONS.
%
%   N is the cutoff frequency, and C must provide the moments of the
%   measure up to order N.
%
%   Supported OPTIONS:
%		- kernel: Fejer(*) | Gaussian | Gaussian-squared | Jackson | K3
%		- grid_size
%		- draw  : 0(*) | 1

%fprintf('Computing convolution proxy ...');

d = numel(n);
N = prod(n+1);

% load options
order  = getoptions(options,'order','colex');
kernel = getoptions(options,'kernel','Fejer');
mgrid  = getoptions(options,'grid_size',101*ones(1,d)); % evaluation grid
M      = prod(mgrid);
draw   = getoptions(options,'draw',0);
nbloc  = getoptions(options,'parallelize',0);

% safeguards
% grid must be thinner than cutoff frequency
gc = ceil((mgrid-1)/2);
fc = n+1;
if sum(gc<fc) % check along each dimension
   cut  = sprintf('%d ', fc);
   grid = sprintf('%d ', gc);
   error('Grid coarser than cutoff frequency: grid=%s, fc=%s', grid, cut);
end

% Jackson order must be even
if strcmp(kernel,'K2')
   if mod(n,2)
      error('Jackson kernel is only defined for even orders');
   end
end

% Fourier matrix
G  = compute_grid(mgrid,'spatial');
%Xf = compute_grid(n,'spectral');
%F  = exp(-2i*pi * reshape(Xf,[N,1,d]) .* reshape(G,[1,M,d]));

% convolution weights
if strcmp(kernel,'Gaussian') | strcmp(kernel,'Gaussian-squared')
   sig = getoptions(options,'gaussian_var',.01/max(n));
   w   = compute_weights(kernel,n,sig);
else
   w   = compute_weights(kernel,n);
end

% fft implementation
%c = M(moms(n,order,1)); % extract moments from moment matrix
c = reshape(c,[2*n+1,1]);
%
S0 = fftshift(w.*ifftshift(c)); % convolution in Fourier
S  = padarray(S0,ceil ((mgrid-1)/2)-n,'pre');
S  = padarray(S, floor((mgrid+1)/2)-n-1,'post');
%TODO: find a more elegant way

% compute proxy (fft)
%P = (2*Fc)^d * real(exp(-2i*pi*Fc*sum(G,d+1)) .* ifftn(S));
shift = reshape(gc,[1,1,d]);
%P = 2^d*prod(gc) * real(exp(-2i*pi * sum(shift.*G, d+1)) .* ifftn(S));
P = prod(mgrid) *   real(exp(-2i*pi * sum(shift.*G, d+1)) .* ifftn(S));

%     %not fft implementation (in 2D, for checking purposes...)
%     [Yf,Xf]   = meshgrid(0:n(1),0:n(2));
%     F     = @(x) exp(-2i*pi*(Xf(:)*x(:,1)' + Yf(:)*x(:,2)'));
%     evalP = @(M) 1/N*real( sum( ((W.*M)*F(G)) .* conj(F(G)), 1) );
%     % compute proxy
%     P2 = evalP(M);
%     P2 = reshape(P2,[mgrid,mgrid]);

% P = P/sum(abs(P(:)));
P = P/sum(P(:)); % WARNING: this result in signed measure for 'Gaussian' convolution

if draw
   clf;
   if d==1
      plot(G,P,'linewidth',3);
      drawnow;

   elseif d==2
      surf(G(:,:,1),G(:,:,2),P,'linestyle','none');
      view(2);
      drawnow;

   else
      warning('Cannot display for dimensions higher than 2');

   end
end

end
    

function w = compute_weights(type,n,varargin)
d = numel(n);

switch type
    case 'F'
       Xf = compute_grid(n,'spectral-sym');
       w  = prod(1 - abs(Xf) ./ reshape(n+1,[ones(1,d) d]), d+1);
		 w  = ifftshift(w);
		 
		 % comparison with auto-correlation computation
% 		 O = ones([n+1,1]); O = padarray(O,n,'post');
% 		 w2 = ifftn(abs(fftn(O)).^2) / prod(n+1);
		 
		 % comparison with fft
%		 L = 256;
%		 [Y,X] = meshgrid((0:2*L-1)'/2/L);
%		 f1 = sin(pi*(n(1)+1)*X).^2 ./ sin(pi*X).^2;
%		 f1(isnan(f1)) = (n(1)+1)^2;
%		 f2 = sin(pi*(n(2)+1)*Y).^2 ./ sin(pi.*Y).^2;
%		 f2(isnan(f2)) = (n(2)+1)^2;
%		 f = f1.*f2;
%		 w_test = real(fftshift(fft2(f)));
%		 w_test = 1/4/L/L * w3(L+1-n(1):L+1+n(1),L+1-n(2):L+1+n(2));
%		 w_test = 1/prod(n+1) * ifftshift(w3);
		 
%        W = ones([2*n+1,1]);
        
    case 'Gaussian'
       Xf  = compute_grid(n,'spectral-sym');
       sig = varargin{1};
       w   = (sqrt(2*pi)*sig)^d * exp(-2*pi^2*sig^2*sum(Xf.^2,d+1));
		 w   = ifftshift(w);
        
%         [V,U] = meshgrid(0:n(1));
%         toepl = @(s) double(U-V)== s; % HACK: does not work if n(1) != n(2)...
%         O = zeros(2*n+1); O(1:n+1,1:n+1) = 1;
%         D = fftshift(ifft2(abs(fft2(O)).^2));
%         G = N*w./D; % HACK: why N?
%         W = zeros(N);
%         for k1=-n(1):n(1)
%             for k2=-n(2):n(2)
%                 th = kron(toepl(k2),toepl(k1));
%                 W = W + G(k1+n(1)+1,k2+n(2)+1)*th;
%             end
%         end
        
    case 'Gaussian-squared'
       Xf  = compute_grid(n,'spectral');
       sig = varargin{1};
       w   = (sqrt(2*pi)*sig)^d * exp(-2*pi^2*sig^2*sum(Xf.^2,d+1)); % w0
       w   = padarray(w,n,'post');
       w   = ifftn(abs(fftn(w)).^2);
        
%         W = N*w0(:).*w0(:)'; %HACK: why N?
        
    case 'K2'
		 m = floor(n/2);
		 O = ones([m+1,1]); O = padarray(O,m,'post'); % outer padding
		 D = fftshift(ifftn(abs(fftn(O)).^2));
		 w = ifftn(abs(fftn(padarray(D,m))).^2); % with centered padding;

		 % comparison with fft approximation
%		 L = 256;
%		 [Y,X] = meshgrid((0:(2*L-1))'/2/L);
%		 %f1 = sin(pi*(2*m(1)-1)*X).^4 ./ sin(pi*X).^4;
%		 f1 = sin(pi*(m(1)+1)*X).^4 ./ sin(pi*X).^4;
%		 f1(isnan(f1)) = (m(1)+1)^4;
%		 f2 = sin(pi*(m(2)+1)*Y).^4 ./ sin(pi.*Y).^4;
%		 f2(isnan(f2)) = (m(2)+1)^4;
%		 f = f1.*f2;
%		 w_test = real(fftshift(fft2(f)));
%		 w_test = 1/4/L/L * w3(L+1-2*m(1):L+1+2*m(1),L+1-2*m(2):L+1+2*m(2));
%		 w_test = ifftshift(w3); % divide by prod(n+1)??
 
%         W = N*D(:).*D(:)';

	 case 'K3'
         m = floor(n/3);

			%TODO: exact formula
%        O1 = ones([m+1,1]); O1 = padarray(O1,m,'post');
%        O2 = fftshift(ifftn(abs(fftn(O1)).^2));
% 			O2 = padarray(O2,m);
%        O3 = ifftn(abs(fftn(O2)).^2);
% 			%O3 = padarray(O3,m,'post');
%        w = ifftn(abs(fftn(O3)).^2);

         % fft approximation
			L = 2048;
			[Y,X] = meshgrid((0:2*L-1)'/2/L);
			f1 = sin(pi*(m(1)+1)*X).^6 ./ sin(pi*X).^6;
			f1(isnan(f1)) = (m(1)+1)^6;
			f2 = sin(pi*(m(2)+1)*Y).^6 ./ sin(pi*Y).^6;
			f2(isnan(f2)) = (m(2)+1)^6;
			f = f1.*f2;
			w = real(fftshift(fft2(f)));
			w = 1/4/L/L * w(L+1-n(1):L+1+n(1),L+1-n(2):L+1+n(2));
			%w = 1/prod(n+1)*ifftshift(w3);
			%w(w<0) = 0;
			w = ifftshift(w);
        
    otherwise
		 
        error(['Unknown kernel: "', type, '"']);
end

%w = ifftshift(w);

%w = ifftshift(w); % put 0 frequency in corner
%fprintf(' Done.\n');
end

