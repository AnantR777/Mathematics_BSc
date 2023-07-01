%% Anantajit Raja
% MATH0033 Numerical Methods Computational homework 3
%
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
% 
f=@(t, y) 1 - y.^2;
tspan = [0 20];
y0 = 0;
N = [19, 21, 40, 80, 160];

% plotting exact solns
for i = 1:5
    t = [0:20/N(i):20];
    y_t = @(t) (exp(2.*t) - 1)/(exp(2.*t) + 1);
end
plot(t, y_t(t));
hold on

% plotting solns of FE methods
for i = 1:5
    diff = zeros(1, N(i));
    [t, u_fe] = feuler(f, tspan, y0, N(i));
    plot(t, u_fe);
for j = 1:length(t)
    diff(j) = abs(u_fe(j) - y_t(t(j)));
end
e_fe(i) = max(diff);
e_feend(i) = max(abs(u_fe(end) - y_t(20)));
end
legend('Exact soln', 'FE N = 19', 'FE N = 21', 'FE N = 40', 'FE N = 80', 'FE N = 160');
title('solns of FE methods vs exact soln')
hold off
e_fe
e_feend
%%
% (b)
% 
f=@(t, y) 1 - y.^2;
tspan = [0 20];
y0 = 0;
N = [19, 21, 40, 80, 160];

% plotting exact solns
for i = 1:5
    t = [0:20/N(i):20];
    y_t = @(t) (exp(2.*t) - 1)/(exp(2.*t) + 1);
end
plot(t, y_t(t));
hold on

% plotting solns of Heun's method
for i = 1:5
    diff = zeros(1, N(i));
    [t, u_heun] = heun(f, tspan, y0, N(i));
    plot(t, u_heun);
for j = 1:length(t)
    diff(j) = abs(u_heun(j) - y_t(t(j)));
end
e_heun(i) = max(diff);
e_heunend(i) = max(abs(u_heun(end) - y_t(20)));
end
legend('Exact soln', 'Heun N = 19', 'Heun N = 21', 'Heun N = 40', 'Heun N = 80','Heun N = 160');
title('solns of Heun methods vs exact soln')
hold off
e_heun
e_heunend
%%
% (c)
% 
N = [19, 21, 40, 80, 160];
loglog(N, e_fe, N, e_heun)
legend('FE method', 'Heun method')
title('e_{fe} and e_{heun} vs N')
% and so it turns out that this agrees with the theory from lectures:
% the Forward Euler method here seems to be first order convergent and
% Heun's method here is second order convergent.
%%
% (d)(i)
% Note from the graph of the FE method we get that the approximations show
% asymptotic behaviour at 1 for N greater than or equal to 21. Therefore the
% FE method well-approximates y(20) for N equal to 21, 40, 80 and 160 (as
% can be seen by the small/zero values of e_feend, but it's not a good
% approximation for N equal to 19. For Heun's method we get that the
% approximations show asymptotic behaviour for N greater than or equal to
% 40. Therefore Heun's method well-approximates y(20) for N equal to 40, 80
% and 160 (which can be seen by the zero values of e_heunend), but it's not a
% good approximation for N equal to 19 or 21.

%%
% (ii)
% Note that the |f_y| is equal to 2y. Note also that y = tanh(t) must also
% be in between 0 and 1 for t > 0 (by examination of the tanh graph). Hence
% |f_y|'s max value is lambda = 2 on the solution trajectory, and so the
% methods should be stable for h < 2/|lambda| = 2/2 = 1, i.e. the critical
% stability value is 1. It turns out that our numerical results agree with
% this because, for the N = 19 case, using the FE method we obtain h =
% 20/19, which is greater than 1 and it's not stable, whereas
% for N greater than 21 we get h < 1 which is stable (as expected).
% Now for Heun's method we get that at N = 19 there's decay to a different
% horizontal asymptotic value and also at N = 21 the error decays to another
% asymptote as t approaches infinity. For N greater than or equal to 40 we
% get stability since h < 1 (as expected).

%%
% (iii)
% For the FE method the lack of stability at N = 19 manifests itself by
% being periodically increasing as t approaches infinity. For Heun's method
% the lack of stability at N = 19 manifests itself by decaying to some
% horizontal asymptote unequal to 1.