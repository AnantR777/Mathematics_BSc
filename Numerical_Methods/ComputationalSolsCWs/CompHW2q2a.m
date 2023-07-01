%% Anantajit Raja
% MATH0033 Numerical Methods Computational homework 2
%
clear all, close all, clc, format long, format compact
%%
% *Exercise 2*
%
% (a)
% (i)
rhoBJ = zeros(4,1); rhoGS=zeros(4,1);
% initialise spectral radii
Nvec = [5,10,20,40,80];
for i = 1:5;
    % iterate for each value in N_values
    N = Nvec(i);
    h = 1/N;
    A= (2/h^2)*diag(ones(N - 1, 1)) - (1/h^2)*diag(ones(N - 2, 1), 1)...
    - (1/h^2)*diag(ones(N - 2, 1), -1);
    D = diag(diag(A));
    L = tril(A)- D;
    U = triu(A) - D;
    % form D, L, U
    B_Jacobi = -(D)^(-1)*(L + U);
    % form B_Jacobi using D, L, U
    B_GS = -(L + D)^(-1)*U;
    % form B_GS using D, L, U
    rhoBJ(i)=max(abs(eig(B_Jacobi)));
    rhoBGS(i)=max(abs(eig(B_GS)));
end
rhoBJ, rhoBGS

% (ii) Yes, both methods should converge for all values of N because we
% have that the spectral radius of the iteration matrices for both methods
% is less than 1. Note that because the spectral radius of the iteration
% matrix for the Gauss-Seidel method is smaller than that for the Jacobi
% method, the Gauss-Seidel method converges faster.

% (iii) Both spectral radii tend to 1 as N increases. In terms of
% performance, both methods will converge slower for larger N.

% (iv)
loglog(Nvec, 1 - rhoBJ)
hold on
loglog(Nvec, 1 - rhoBGS)
loglog(Nvec, Nvec.^(-2))

xlabel('N')
ylabel('1-\rho(B)')
legend('Jacobi Method','Gauss-Seidel Method', 'N^{-2}', 'Location','NorthEast')

% (v)
% From the plot we have $log(1-\rho (B))$ vs $log(N)$ is linear and so we
% have $1-\rho (B)) \approx CN^\alpha$, i.e. $\rho (B) \approx 1 -
% CN^\alpha$. Now we can find $\alpha$ and C using (rearranging the $\rho$
% formula):

J_grad = (log(1 - rhoBJ(5)) - log(1 - rhoBJ(3)))/(log(Nvec(5)) - log(Nvec(3)))
GS_grad = (log(1 - rhoBGS(5)) - log(1 - rhoBGS(3)))/(log(Nvec(5)) - log(Nvec(3)))
C_J = Nvec(5)^-J_grad*(1 - rhoBJ(5))
C_GS = Nvec(5)^-GS_grad*(1 - rhoBGS(5))

% giving $\alpha \approx -2$; for the Jacobi method $C \approx 4.9$ but
% for the Gauss-Seidel method $C \approx 9.6$.

%(vi)
% We can deduce err(k) = \beta (1-CN^\alpha)^kerr(0) for some \beta by
% the proportionality statement given in the question. Hence err(k) = \beta
% (1-CN^\alpha)^kerr(0) and so putting h = 1/N we get err(k) = \beta
% (1-Ch^2)^kerr(0). Using the fact from notes that err(k)/err(0) = tol gives 
% (1-Ch^2)^k = \frac{tol}{\beta} and so k*log(1-Ch^2) =
% log(\frac{tol}{\beta}). Using small order approximation for h we get
% -Ch^2k = log(\frac{tol}{beta}) and therefore k =
% [log(\frac{\beta}{tol})]N^2. We must have that k is O(N^2) so that the
% coefficient of N^2 in the k-equation is constant (and hence the error
% too). Hence we get that the total cost for these methods is O(N^3), that
% is the same as direct methods (we combine the fact for k with 
% the fact that cost is O(n) for a tridiagonal matrix A being multiplied
% by a vector).