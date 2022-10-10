function Tv = Tprod(c,v,varargin)
% TPROD perform Toeplitz tensor multiplication

%TODO: add flags on sizes!
if ndims(c) ~= ndims(v) %TODO: doesn t work for 2-d vs 1-d...
	% check for obvious incorrect formatting of inputs
	error('Dimension mismatch between tensors');
end

crop_size = cellfun(@(i) (1:i), num2cell(size(v)), 'uniformOutput', false);

if nargin > 2
	tflag = varargin{1};
else
	Tv = ifftn(fftn(c) .* fftn(v,size(c))); % Toeplitz-vector = convolution
	Tv = Tv(crop_size{:}); % crop
	Tv = Tv(:);
	return;
end

% specifically for use in svds
if strcmp(tflag,'notransp')
	Tv = ifftn(fftn(c) .* fftn(v,size(c)));
else
	Tv = ifftn(fftn(c) .* fftn(v,size(c))); %TODO: handle non Hermitian cases
end
Tv = Tv(crop_size{:}); % crop
Tv = Tv(:);
