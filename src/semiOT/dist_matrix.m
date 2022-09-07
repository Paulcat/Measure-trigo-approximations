function C = dist_matrix(x,y,p)
%COST_MATRIX Cost matrix for optimal transport
%   C = COST_MATRIX(X,Y,P) is the matrix with entries Cij =  ||Xi-Yj||_P

%fprintf('\t Computing distance matrix ...');

d = size(x,2); % dimension

% safety nets
if d~=size(y,2)
	error( ['Dimension mismatch between the two measures: ', int2str(d), ...
		' vs ', int2str(size(y,2))] );
end

m = size(x,1);
n = size(y,1);

x = reshape(x,[m,1,d]);
y = reshape(y,[1,n,d]);

%if ~nbloc
D = min(abs(x-y),1-abs(x-y)); % wrap-around component-wise distance
C = sum(D.^p,3).^(1/p);
%else
%	B = bloc_divide(m,nbloc);
%	C = [];
%	parfor i=1:length(B)-1
%		bi = B(i):(B(i+1)-1);
%		xi = x(bi,:,:);
		%
%		Di = min(abs(xi-y),1-abs(xi-y));
%		C  = [C; sum(Di.^p,3).^(1/p)];
%	end
%end

%fprintf(' Done.\n');

end

