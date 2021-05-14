%%
%MATH0033 Numerical Methods Computational Homework 1
%Shiqi Su 
%% 
% Set up  
clear all, close all,clc
format long, format compact
fs=16;
set(groot,'defaulttextfontsize',fs);
set(groot,'defaultaxesfontsize',fs);
set(groot,'defaultLineLineWidth',2)      
set(groot,'defaultContourLineWidth',2) 
set(0,'DefaultLegendAutoUpdate','off') 
%% 
% Excersise_1_(a)
% Bisection Method choose interval [-pi/2,3]
nmax=30;
tol=1e-10;

f=@(x)x/2-sin(x)+pi/6-sqrt(3)/2;
x=linspace(-pi/2,3,100);
figure
plot(x,f(x))
grid on
xlabel('x')
ylabel('f(x)')
title('bisection method')
hold on

[zero,res,niter,itersb]=bisection(f,-pi/2,3,tol,nmax)
for i=1:10
    % Plot the first 10 iterates from bisection
    scatter(itersb(i,1),f(itersb(i,1)),'r')
    % pausing after each one
    pause                                     
end
%%
% To conclude, one of the roots approximated by bisection method is about 2.2460'
%%
% Excersise_1_(b) Newton Method
tol=1e-10;
df=@(x)1/2-cos(x);
f=@(x)x/2-sin(x)+pi/6-sqrt(3)/2;
tol=1e-10;
nmax=30;
x_alpha=pi;
x_beta=-pi/2;
[zero,res,niter,itersn]=newton(f,df,x_alpha,tol,nmax)
[zero,res,niter,itersn]=newton(f,df,x_beta,tol,nmax)
%%
% The solution of f(x) is \alpha=2.24601 and \beta=-1.04720. And the iteration
% for \alpha is 5, for \beta is 27. The difference between to roots is because
% of difference convergence rate. The function is twice continuously
% differentiable at f(\alpha),ie.$df/dx(beta)0$ which implies quadratic convergence. 
% But at point beta,it appears df/dx(beta)~=0 which means liner convergence.
%% 
% Excersise_1_(c) Modified Newton Method
f=@(x)x/2-sin(x)+pi/6-sqrt(3)/2;
df=@(x)1/2-cos(x);
phi=@(x)x-2*f(x)/df(x);
x0=-pi/2;
[fixp,res,niter,itersfp1]=fixpoint(phi,x0,tol,nmax)
%%
%the number of iteration is 4 
%% 
% Excersise_2_(a)
%(a) define function and its derivative
f=@(x)x+exp(-20.*x.^2).*cos(x);
df=@(x)1+(-exp(-20*x^2)*sin(x)-cos(x)*exp(-20*x^2)*40*x);
%the range of x
x=linspace(-1,1,100);
%plot f(x)
figure
plot(x,f(x),'b','LineWidth',2)
grid on
xlabel('x');
ylabel('f(x)');
title('Exact Function')
hold on

% Apply Netwon Method to get first 10 iterations
tol=1e-10;
nmax=30;
x0=0;
[zero,res,niter,itersn]=newton(f,df,x0,tol,nmax)
disp('value of 10 first iteration')
disp(itersn(1:10))
%% 
% what happens when k tends to infinity
f=@(x)x+exp(-20.*x.^2).*cos(x);
df=@(x)1+(-exp(-20*x^2)*sin(x)-cos(x)*exp(-20*x^2)*40*x);

x=linspace(-1,1,100);
% Apply Netwon Method to get first 10 iterations
tol=1e-10;
nmax=500;
[zero,res,niter,itersn]=newton(f,df,x0,tol,nmax);  
%%
% It appears with oscillation and doesn't converge to roots. 
% Because its convergence theory is for "local" convergence with initial guess, $x_0\in[\alpha-\delta,\alpha+\delta]$
% where $\exists \delta$ is sufficently small.
% Far away from the root the method can have highly nontrivial dynamics.
%% 
% Excersise_2_(b)
% using bisection method to find a better starting point
f=@(x)x+exp(-20.*x.^2).*cos(x);
df=@(x)1+(-exp(-20*x^2)*sin(x)-cos(x)*exp(-20*x^2)*40*x);

nmax=5;
tol=1e-10;
[zero,res,niter,itersb]=bisection(f,-1,1,tol,nmax)

x0=zero;
[zero,res,niter,itersn]=newton(f,df,x0,tol,nmax)
%% 
% In conclusion, the root is -0.257298 with 4 iterations.