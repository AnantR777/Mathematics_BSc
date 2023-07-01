%% Anantajit Raja
% MATH0033 Numerical Methods Computational homework 1
%
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
f = @(x)x/2 - sin(x)+ pi/6 - sqrt(3)/2;
[zero, res, niter] = bisection(f, 1, 3, 1e-10, 100)
%%
% (b)
f = @(x)x/2 - sin(x)+ pi/6 - sqrt(3)/2;
df = @(x)1/2-cos(x);
% for x_alpha = pi
[zero, res, niter] = newton(f,df,pi,10e-10,100)
% for x_beta = -pi/2
[zero, res, niter] = newton(f,df,-pi/2,10e-10,100)
% Note that the more negative root requires more iterations (27) because there
% we have a local maximum and hence linear convergence to that root,
% compared to the more positive root to which there's quadratic
% convergence and hence less iterations (5)
%%
% (c)
f = @(x)x/2 - sin(x)+ pi/6 - sqrt(3)/2;
df = @(x)1/2-cos(x);
phi=@(x)x-2*f(x)./df(x);
% for x_beta = -pi/2
[p,res,niter]=fixpoint(phi,-pi/2,1e-10,100)
