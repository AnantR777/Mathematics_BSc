%% Setup
setup;
clf

A=[2,3;1,4];    %Define the matrix A
b= [2;1]; % Right-hand side vector b
x= A\b;         %Find the "exact solution" (to machine precision) using direct solver
x0 = zeros(2,1); %Set x0 zero initial guess
%x0 = rand(2,1); %Set x0 random initial guess
Nmax = 100; % Set maximum number of iterations
p = 2; % Choose the p-norm used to measure the error. Eg. p=2, p=1 or p=Inf


D=diag(diag(A)); %Extract the diagonal matrix of A

itersJ=zeros(2,Nmax);
itersJ(:,1)=x0;

for k=1:Nmax-1
    %Jacobi  x_(k+1) = x_k + D^-1 (b-A*x_k)
    itersJ(:,k+1) = itersJ(:,k) + inv(D)*(b-A*itersJ(:,k)); % Better implementations avoid using inv(D)
end

errsJ=vecnorm(x-itersJ,p);

semilogy(1:Nmax,errsJ,'*-r');

legend('Jacobi');


