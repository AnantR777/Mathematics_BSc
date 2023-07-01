%% Anantajit Raja
% MATH0033 Numerical Methods Computational homework 2
%
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
% 
% Note since the main diagonal consists of just 1s, we can find the row
% with the strictest condition such that the entry in the diagonal of the
% matrix for that row is bigger than the absolute sum of all the other
% entries in the same row. This 'strictest' condition can be found for the
% middle n - 2 rows of the matrix A_eps:
% 1 > 2(eps + eps ^2)
% which implies 0 < eps < (\sqrt{3} - 1)/2 = 0.366 (3sf)
%%
% (b)
% 
[A,b] = matrix(5,0.3);
% n = 5 and eps = 0.3
x_0 = zeros(5,1);
% n = 5 column consisting of 1s
[x_vec, JacobiIter]= itermeth(A, b, x_0, 1e3, 1e-10, 'J')
% nmax = 1000, tol = 10^-10
[x_vec, GSIter]= itermeth(A, b, x_0, 1e3, 1e-10, 'G')

% i.e. the Jacobi method converged in 50 iterations and the Gauss-Seidel
% method converged in 14 iterations.
%%
% (c)
% 
num = 101;
% number of epsilons to plot
spectralJacobi = zeros(num, 1); spectralGS = zeros(num, 1);
% initialises matrices
eps = linspace(0, 1, num);
% need to plot for eps = 0, 0.01, ... , 1
for i = 1:num;
    A = matrix(5,eps(i));
    % constructs A and b
    D = diag(diag(A));
    L = tril(A)- D;
    U = triu(A) - D;
    % form D, L, U
    B_Jacobi = -(D)^(-1)*(L + U);
    % form B_Jacobi using D, L, U
    B_GS = -(L + D)^(-1)*U;
    % form B_GS using D, L, U
    spectralJacobi(i) = max(abs(eig(B_Jacobi)));
    % finds spectral radius of B_Jacobi
    spectralGS(i) = max(abs(eig(B_GS)));
    % finds spectral radius of B_GS
end
plot(eps,spectralJacobi,'b-')
% plots rho vs eps for Jacobi
hold on
plot(eps,spectralGS,'m--')
% plots rho vs eps for GS
xlabel('\epsilon'), ylabel('Spectral Radius \rho')
legend('Jacobi Method','Gauss-Seidel Method','Location','NorthWest')
hold off

% (i) Note a method converges iff the spectral radius of the iteration matrix 
% is less than 1. Hence the Jacobi method here converges for values of 
% epsilon less than about 0.4 and the Gauss-Seidel method here converges
% for values of epsilon less than about 0.7. The answer obtained here is
% about double the upper bound for the value of epsilon obtained in part (a).

% (ii) The Gauss-Seidel method converges faster because the
% spectral radius of that iteration matrix is always lower than that of the
% spectral radius of the iteration matrix for Jacobi method.

% (iii) The Gauss-Seidel method will converge for $\epsilon = 0.5% since
% the spectral radius of B_GS is less than 1, however the Jacobi
% method will not converge for this %\epsilon$ because the spectral radius
% is greater than 1 here. These can be found by:

[A,b] = matrix(5, 0.5);
% n = 5 and eps = 0.5
x_0 = zeros(5, 1);
% n = 5 column consisting of 1s
[x_vec, JacobiIter] = itermeth(A,b,x_0,1e3,1e-10,'J')
% nmax = 1000, tol = 10^-10
[x_vec, GSIter] = itermeth(A,b,x_0,1e3,1e-10,'G')