%(b), (c)
Nvec = [5, 10, 20, 40, 80];
hvect = 1./Nvec;
% h = 1/N, elementwise
error_vect1 = zeros(5, 1);
% same length as Nvec, initialise error vector
error_vect2 = zeros(5, 1);
close all
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
error_vect1
error_vect2
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