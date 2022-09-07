function [X,A,P] = multiscale_decomposition(x0,a0,L,S,p,w)
%MULTISCALE_DECOMPOSITION Compute multiscale decomposition of a measure
%	[Y,B,P] = MULTISCALE_DECOMPOSITION(X,A,L,S,P,W) returns a (L+1)-level
%	decomposition Y = (Y0, ..., YL) of X, such that
%       - Y0 = X
%       - Y(i) contains S(i) points
%       - Y(i+1) has a low Wasserstein distance to Y(i)
%   
%   This is performed using Lloyd's algorithm. P(i) contains the transport
%   plan between Y(i) and Y(i+1), i.e. ...
%
%   Q specifies the norm and W the Wasserstein distance used.

%disp('Computing quantization...');

% initialize;
X = cell(1,L+1); X{1} = x0;
A = cell(1,L+1); A{1} = a0;
P = cell(1,L); % if L=0, P is empty


x 	 = x0;
a 	 = a0;
tol = 1e-2;


for l=1:L
	% random subset of points
	I = randperm(size(x,1),S(l+1));
	y = x(I,:);
	
	[x,a,J] = lloyd(y,x,a,tol,p,w,0);
		
	% store
	X{l+1}  = x;
	A{l+1}  = a;
	P{l} 	  = J;
end

%disp('...end');
end
