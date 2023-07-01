phi1=@(x)1-sin(2*x);
phi2=@(x).5*asin(1-x);
[fixp,res,niter]=fixpoint(phi1,0.7,1e-8,1000)
[fixp,res,niter]=fixpoint(phi2,0.7,1e-8,1000)