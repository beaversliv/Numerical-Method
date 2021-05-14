% D. Hewett MATH0033 Numerical Methods 
% demo_SIR.m
% Demonstration of the basic SIR model for infectious diseases. 
% The system is
% dS/dt = -beta*S*I
% dI/dt = beta*S*I - gamma*I
% dR/dt = gamma*I
% Here S, I and R represent the proportions of the population that are
% Susceptible, Infected and Removed (i.e. recovered or died). 
% 1/gamma = typical time from infection until removal (recovery or death)
% 1/beta = typical time between contacts
%
% This is a 3-dimensional system of the form 
%      y'(t)=f(t,y(t))
% where y=(S;I;R) (column vector) 
% and f=(-beta*S*I;beta*S*I-gamma*I;gamma*I) (column vector)
%
% In this demo code we show how a one-step method (Heun's method) can be
% used to solve an IVP for this system, given choices of beta, gamma, and 
% the initial values of S, I and R.
%
% Here we actually specify gamma and the ratio R0:=beta/gamma, and then
% compute beta as beta=R0*gamma.
% 
% Try running the code for different choices of R0 (either <1 or >1) and
% different initial data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose parameters in the model (ask your favourite epidemiologist)
gamma=1/7; % 1/gamma = typical time from infection until removal (recovery or death)
R0=2; % R_0= beta/gamma. 
y0=[0.9;0.1;0]; % Initial values of S, I and R
tmax=60;    % Length of the time interval over which to solve
M=100;     % Number of mesh elements in the discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solution of the IVP
beta=R0*gamma;  % 1/beta = typical time between contacts
f=@(t,y)[-beta*y(1)*y(2); beta*y(1)*y(2)-gamma*y(2); gamma*y(2)];
t = linspace(0,tmax,M+1);
h = tmax/M;
u = zeros(3,M+1);
u(:,1) = y0;
for n = 1:M
    % Heun
    u(:,n+1)=u(:,n)+h/2*(f(t(n),u(:,n))+f(t(n+1),u(:,n)+h*f(t(n),u(:,n))));
    % FE (in case you want to play with that)
%    u(:,n+1)=u(:,n)+h*(f(t(n),u(:,n)));
end
figure
plot(t,u(1,:),'kx-',t,u(2,:),'rx-',t,u(3,:),'gx-')
legend('S','I','R')
xlabel('time (days)')
ylabel('proportion of population')
title(['SIR model with \gamma=' num2str(gamma) ' and R_0=' num2str(R0)])