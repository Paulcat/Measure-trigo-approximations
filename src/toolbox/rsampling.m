function x = rsampling(f,M,s,a,b)
%RSAMPLE Rejection sampling with uniform bound

y = unifrnd(a,b,10*s,1);
u = rand(10*s,1);

x = y(u < (f(y)/M));
x = x(1:s);

