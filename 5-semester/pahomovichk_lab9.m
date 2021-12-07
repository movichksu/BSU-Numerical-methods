clear;
clc;

f=@(x) x.^2 - 5*x + 6;
tau = 1/9;
eps = 10**(-11);
fi=@(x) x+tau*f(x); 
X = [1.5:0.1:2.5];

figure
plot(X, X, X, fi(X))

errs = [];
x = X;

x(1) = 1.5;
k=1;
Err=1;
errs(1) = Err;
while Err > eps  && k < 100
  x(k+1) = fi(x(k));
  Err = abs(f(x(k+1)));
  errs(1,k) = x(k+1);
  errs(2,k) = Err;
  k = k + 1;
endwhile

sprintf('answer x=%0.11f',x(k-1))
semilogyerr (errs(1, :), errs(2, :));