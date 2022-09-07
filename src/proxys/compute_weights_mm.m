function W = compute_weights_mm(type,n,varargin)
d = numel(n);

switch type
    case 'Fejer'
       Xf = compute_grid(n,'spectral-sym');
       w  = prod(1 - abs(Xf) ./ reshape(n+1,[ones(1,d) d]), d+1);
		 %w  = ifftshift(w);
       %W = ones([2*n+1,1]);
		 W = ones(prod(n+1));
        
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
        
        W = prod(n+1)*w(:).*w(:)'; %HACK: why N?
        
    case 'Jackson'
        m = floor(n/2);
        O = ones([m+1,1]); O = padarray(O,m,'post'); % outer padding
        D = fftshift(ifftn(abs(fftn(O)).^2));
        %w = ifftn(abs(fftn(padarray(D,m))).^2); % with centered padding;
		  
        
        %W = prod(n+1)*D(:).*D(:)';
		  W = D(:).*D(:)';
        
    otherwise
        error(['Unknown kernel: "', type, '"']);
end

%w = ifftshift(w);

%w = ifftshift(w); % put 0 frequency in corner
%fprintf(' Done.\n');
end