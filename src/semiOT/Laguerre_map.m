function L = Laguerre_map(D,f)
%LAGUERRE_MAP Compute Laguerre cells
%   L = LAGUERRE_MAP(D,F,DIM) computes the Laguerre cells associated to the
%   weights F, for the metric D, over some grid
%
%   D is a matrix such that D(i,j) = dist(i,j), for i,j in some grid I 
%   L is a vector such that L(i) = k s.t. i is in Laguerre_k(F), for i in I

[~,L] = min(D-f(:)',[],2);
%L = reshape(L,[m,m]);

end

