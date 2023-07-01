f = @(x)sin(2*x)-1+x;
x=[-3:0.1:3];
plot(x,f(x))
grid on

[zero,res,niter]=bisection(f,-1,1,1e-8,1000)

% execute publish command from command window, not in this script