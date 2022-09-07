function X = repsym_curve(x,w)
%WARNING: Hardcoded... 

m = length(x);
if mod(m,4)
	error('sampling cannot be divided in 4 parts');
end

s = m/4;
y = zeros(s,4);
for i=1:4
	I = (i-1)*s+1 : i*s;
	y(:,i) = x(I);
end

X = [y(:,1), w(y(:,1)); w(y(:,2)), -y(:,2); -y(:,3), -w(y(:,3)); -w(y(:,4)), y(:,4)];
end
