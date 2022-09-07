function [x,a,J] = lloyd(x0,target,b,tol,p,q,display)
% LLOYD
%	x0 initial centroids, target is measure wrt which Voronoi cells are
%	computed

%TODO: guarantee that the output J match all points of x (no empty Voronoi cells...)

n = size(target,1);
d = size(target,2);

% initialize
m = size(x0,1);
x = x0;
%a = a0;
a = ones(m,1)/m;

maxit = 100;

it		= 0;
crit 	= 1;
while crit >= tol & it<maxit
	% compute Voronoi cells
	D = dist_matrix(target,x,p).^q;
	[~,J] = min(D,[],2);
	
	% barycenter of cells
	A = (J==reshape(1:m,[1 1 m])) .* b;
	B = A .* reshape(target,[n 1 1 d]);
	
	% remove points with empty Voronoi cells?TODO
	E = sum(A,[1,2])==0; %TODO: squeeze or flatten for improve efficiency?
	A(:,:,E) = [];
	B(:,:,E,:) = [];
	
	% update centroids to barycenters
	a_old = a;
	a 		= sum(A,1); % this is normalized
	x_old = x;
	x 		= sum(B,1) ./ a;
	%
	%x 		= squeeze(x);
	x 		= reshape(x,[],d);
	a 		= squeeze(a);
	%
	m_old = m;
	m 		= size(x,1);

	% update criterion: centroid displacement
	crit 	= norm(x_old(~E,:) - x,'fro')/norm(x,'fro');
	it 	= it + 1;

	if display
		clf, hold on;
		j = 1;
		ms1 = min(5/m/max(a)*1e4,40/max(a));
		ms2 = min(5/n/max(b)*1e6,1/max(b));
		for i=1:m_old
			col = rand(1,3);
			I = find(J==i);
			if isempty(I) % Voronoi cell is empty, point removed
				scatter(x_old(i,1),x_old(i,2),ms1*a_old(i),'k.');
			else
				scatter(target(I,1),target(I,2),ms2*b(I),col,'filled');
				scatter(x_old(i,1),x_old(i,2),ms1*a_old(i),'MarkerFaceColor',col,...
					'MarkerEdgeColor','k');
				scatter(x(j,1),x(j,2),5*ms1*a(j),'pentagram','MarkerFaceColor',...
					col,'MarkerEdgeColor','k');
				j = j+1;
			end
	    end
		 pause(.5);
	end


end

%TODO: clarifier la parallelisation
%TODO: etre sur d'avoir bien compris l'algorithme de Lloyd --> ok

end
