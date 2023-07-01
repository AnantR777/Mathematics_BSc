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

%(b), (c)
Nvec = [5, 10, 20, 40, 80];
hvect = 1./Nvec;
% h = 1/N, elementwise
error_vect1 = zeros(5, 1);
% same length as Nvec, initialise error vector
error_vect2 = zeros(5, 1);
close all
% closes previous plots
figure(1), hold on, xlabel('x'), ylabel('y')
figure(2), hold on, xlabel('x'), ylabel('y')
for i = 1:5
    % number of iterations equal to length of Nvec
    N = Nvec(i);
    h = hvect(i);
    A = (2/h^2)*diag(ones(N - 1, 1)) - (1/h^2)*diag(ones(N - 2, 1), 1)...
    - (1/h^2)*diag(ones(N - 2, 1),-1);
    % same definition as before
    x_vec = transpose(linspace(h, 1 - h, N - 1));
    % from 1/N to (N - 1)/N
    b1_vec = transpose(sin(pi*h*(1:N-1)));
    % given in Q, need column vector
    b2_vec = ones(N - 1, 1);
    % column vector of N - 1 1s
    u1_vec = itermeth(A, b1_vec, zeros(N - 1, 1), 1e5, 1e-10, 'G');
    % column of zeros is 
    u2_vec = itermeth(A, b2_vec, zeros(N - 1, 1), 1e5, 1e-10, 'G');
    % our initial vec
    y1_vec = pi^-2*sin(pi*x_vec);
    y2_vec = 1/2*x_vec.*(1 - x_vec);
    % element-wise multiplication
    figure(1), plot([0; x_vec; 1],[0; u1_vec; 0])
    % plots approximations against x
    figure(2), plot([0; x_vec; 1],[0; u2_vec; 0])
    n = 1;
    possible_e1_i = [];
    % take max over this list at the end
    possible_e2_i = [];
    while n <= N - 1
        % appends to lists above
        possible_e1_i(n) = abs(u1_vec(n) - y1_vec(n));
        possible_e2_i(n) = abs(u2_vec(n) - y2_vec(n));
        if n == N - 1
            break
        end
        n = n + 1;
    end
    error_vect1(i) = max(possible_e1_i);
    % finds max over lists above to find entries in error vecs
    error_vect2(i) = max(possible_e2_i);
end
figure(1), title('f = sin(\pix)'), plot([0; x_vec; 1], [0; y1_vec; 0]),
legend('N = 5','N = 10','N = 20','N = 40','N = 80','y(x)')
figure(2), title('f = 1'), plot([0; x_vec; 1], [0; y2_vec; 0]),
legend('N = 5','N = 10','N = 20','N = 40','N = 80','y(x)')
% exact y only plotted once so not in loop
figure
loglog(hvect, error_vect1, 'r-', hvect, error_vect2)
slope1 = (log(error_vect1(5))-log(error_vect1(4)))/(log(hvect(5))-log(hvect(4)))
slope2 = (log(error_vect2(5))-log(error_vect2(4)))/(log(hvect(5))-log(hvect(4)))
C1 = hvect(5)^-slope1*error_vect1(5)
C2 = hvect(5)^-slope2*error_vect2(5)
xlabel('h')
ylabel('error')
hold on
loglog(hvect, C1*hvect.^slope1, 'b--')
loglog(hvect, C2*hvect.^slope2, 'b--')
legend('f = sin(\pix)', 'f = 1', 'linear fit-sin', 'linear fit-1', 'Location', 'NorthWest')

% Note from the plot of error vs h, the error for the sin-graph clearly
% coincides with a linear fit for the graph and so, since the gradient of
% the slope for this is 2, we have p = 2 and so convergence is of $O(h^2)$
% with const. C approximately 0.08, agreeing with the theoretical estimate.

% The thing that's surprising about the error for the 1-graph is that it
% looks like the error is diverging as the error dips below the straight
% line, but since y''''(x) = 0 we end up with C = 0 and so 0 error, i.e.
% an exact solution. Hence the given plot for the error for the 1-graph
% contradicts the true nature of the error, probably because of imprecision
% in coding and solving this linear system.