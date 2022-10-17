function [x,a,mvec] = instantiate_measure(type,s,d,varargin)

normalize1 = @(x) x/norm(x,1);

if nargin < 4
    sampling = 'uniform';
else
    sampling = varargin{1};
end

switch type
	case 'd'
		x = .1 + .8*rand(s,d);
		a = normalize1(.1+.6*rand(s,1));

		mvec = @(kv) exp(-2i*pi*sum(reshape(kv,[],1,d).*reshape(x,[1,s,d]),3))*a;

	case 'line'
        error('TODO');
		if d~=2
			error('my line only works in 2 dimensions');
		end
		q1 = 2; q2 = 3; l = q1/q2; % rational coefficient
		b  = .3; 
		x  = linspace(0,q2,s)';
		y  = l*x + b;
		x  = [mod(x,1), mod(y,1)];
		a  = normalize1(ones(s,1));

		mvec = @(kv) exp(-2i*pi*kv(:,2)*b) .* ((kv(:,1)*q2 + kv(:,2)*q1)==0);

	case 'c'
		if d~=2
			error('my circle only works in 2 dimensions');
		end
		c0 = [.5 .5];
		r0 = .3;

		switch sampling
			case 'uniform'
				theta = 2*pi*rand(s,1);
			case 'regular'
				theta = 2*pi*linspace(0,1,s)'; % list of angles
			otherwise
				error('Sampling can be either "uniform" or "regular"');
		end
		x 		= c0 + r0 * [cos(theta), sin(theta)];
		a		= normalize1(ones(s,1));

		mvec = @(kv) exp(-2i*pi*sum(kv.*c0,d)) .* besselj(0,-2*pi*r0*vecnorm(kv,2,2));

	case 'z'
		if d~=2
			error('my curve only works in 2 dimensions');
		end

		if mod(s,4)
			error('please choose a number of points divisible in 4 parts...');
		end

		c0 = [.5 .5];
		s1 = s/4;

		% parameterization
		w = @(t) acos(-1 + 2*5/8 ./ (cos(2*pi*t) + 1))/2/pi;

		% rejection sampling
		f = @(t) ...
			1/4 * sqrt( ...
						abs( ...
							exp(1i*acos(5./4./(cos(2*pi*t)+1)-1)/2/pi) .* ...
							sin(2*pi*t) ./ ...
							( sqrt(-1/16*(5./(cos(2*pi*t)+1) -4).^2 + 1) .* ...
								(cos(2*pi*t)+1).^2 )...
						).^2 + 16);
		t0 = -acos(-1 + 1/2*sqrt(8*5/8))/2/pi;
		t1 =  acos(-1 + 1/2*sqrt(8*5/8))/2/pi;

		t = zeros(s,1);
		switch sampling
			case 'uniform'
				for i=1:4
					t(((i-1)*s1+1):(i*s1)) = rsampling(f,1.02,s1,t0,t1); %TODO: softcode the upper bound
				end

			case 'regular'
				for i=1:4
					%step = (t1-t0)/s1;
					%t(((i-1)*s1+1):(i*s1)) = t0 + (0:s1-1)'*step;
                    
					% with pre-computed points
					filename = ['points_curve_M',int2str(s),'_quarter.txt'];
					t(((i-1)*s1+1):(i*s1)) = load(filename)/2/pi;

               
					% with Newton iterations
					%step = integral(f,t0,t1)/s1;
					%al   = (0:s1-1)'*step;
					%ti = apply_Newton(f,al,t0,t1);
					%t(((i-1)*s1+1):(i*s1)) = ti;
				end

			otherwise
				error('Sampling can be either "uniform" or "regular"');
		end
		x = repsym_curve(t,w);
		x = c0 + x;
		a = normalize1(ones(s,1));

		mo   = load('moments_curve_n1000.txt');
		mvec = @(kv) mo(sub2ind([1001 1001],abs(kv(:,1))+1,abs(kv(:,2))+1)); %TODO: really symmetric?

	case 'continuous'
		deg  = 40*ones(1,d);
		Xdeg = compute_grid(deg,'spectral-sym');
		Xdeg = reshape(Xdeg,[prod(2*deg+1),d]);
		h 	  = 1./(vecnorm(Xdeg,2,2).^2 + 1);

		tgrid = [100 100];
		s 		= prod(tgrid);
		x 		= compute_grid(tgrid,'spatial');
		x 		= reshape(x,[prod(tgrid),d]);

      B = bloc_divide(prod(tgrid),100);
      a = [];
      for i=1:length(B)-1
			bi = (B(i):B(i+1)-1)';
         ai = exp(+2i*pi * (x(bi,1)*Xdeg(:,1)' + x(bi,2)*Xdeg(:,2)')) * h(:);
         a  = [a;ai];
      end
		a = real(normalize1(a));

		mvec = @(kv) 1./(vecnorm(kv,2,2).^2+1);

	otherwise
		error('I do not know this type of measure, sorry...');
end
