%% Setup
setup;
clf

A=[2,3;1,4];    %Define the matrix A
b= [pi;-sqrt(2)]; % Right-hand side vector b
x= A\b;         %Find the "exact solution" (to machine precision) using direct solver
x0 = zeros(2,1); %Set x0 zero initial guess
%x0 = rand(2,1); %Set x0 random initial guess
Nmax = 100; % Set maximum number of iterations
p = 2; % Choose the p-norm used to measure the error. Eg. p=2, p=1 or p=Inf

%% Richardson Method

alpha_1 = 1/5; 
alpha_2 = 1/4;
alpha_3 = 1/3; % This is the optimal value of the parameter in Richardson method

itersb_1=zeros(2,Nmax);     %Create arrays for all the iterates for the Richardson method
itersb_2=zeros(2,Nmax);
itersb_3=zeros(2,Nmax);

itersb_1(:,1)=x0;       %If using random initial guess, set first iterate to x0.
itersb_2(:,1)=x0;
itersb_3(:,1)=x0;

for k=1:Nmax-1
    
    % Richardson iteration
    % x_(k+1) = x_k + alpha *(b-A*x_k)

    itersb_1(:,k+1) = itersb_1(:,k) + alpha_1*(b-A*itersb_1(:,k));
    itersb_2(:,k+1) = itersb_2(:,k) + alpha_2*(b-A*itersb_2(:,k));
    itersb_3(:,k+1) = itersb_3(:,k) + alpha_3*(b-A*itersb_3(:,k));

    % In Matlab: M(:,k) is the k-th column of the matrix M, so here
    % itersb_1(:,k) denotes the (k-1)-th iteration

end

% Compute the 2 norm errors for each iteration

err_1 = vecnorm(x-itersb_1,p);
err_2 = vecnorm(x-itersb_2,p);
err_3 = vecnorm(x-itersb_3,p);

semilogy(1:Nmax,err_1,'-r',1:Nmax,err_2,'-g',1:Nmax,err_3,'-b');

hold on

%% Jacobi and Gauss-Seidel

D=diag(diag(A)); %Extract the diagonal matrix of A
LD = tril(A);   %Extract the diagonal and lower triangular part of A: LD=L+D

itersJ=zeros(2,Nmax);
itersJ(:,1)=x0;
itersGS=zeros(2,Nmax);
itersGS(:,1)=x0;

for k=1:Nmax-1
    %Jacobi  x_(k+1) = x_k + D^-1 (b-A*x_k)
    itersJ(:,k+1) = itersJ(:,k) + inv(D)*(b-A*itersJ(:,k)); % Better implementations avoid using inv(D)
    %Gauss-Seidel x_(k+1)= x_k + (L+D)^-1 (b - A*x_k)
    itersGS(:,k+1) = itersGS(:,k) + inv(LD)*(b-A*itersGS(:,k)); % Better implementations avoid using inv(LD)

end

errsJ=vecnorm(x-itersJ,p);
errsGS=vecnorm(x-itersGS,p);

semilogy(1:Nmax,errsJ,'--r',1:Nmax,errsGS,':c');

str1 = 'Richardson 1';
str2 = 'Richardson 2';
str3 =  'Richardson 3';
str4 = 'Jacobi';
str5 = 'GS';

legend(str1,str2,str3,str4,str5);
hold off

