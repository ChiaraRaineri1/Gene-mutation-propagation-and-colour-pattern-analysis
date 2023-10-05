%% PROJECT: Fisher equation %%

clc
clear all
close all

%% plot u0

a = -4;
b = 4;
t0 = 0;
T = 10;

u0 = @(x) 0.5 * (cos(pi*(x./2 + 3/4))).^4 .* (x<=1.5) .* (x>=-0.5);

x_disp = linspace(a, b, 1000);

figure(1)
plot(x_disp, u0(x_disp), 'r', 'linewidth', 2);
xlabel('Population');
ylabel('u0');
title('Initial profile of the spreading of the mutation');
grid on; box on;

% max probability (0.5) in 0.5
% start: very few elements of the pop manifest the mutation

%% approximation

% method = UNCONDITIONALLY STABLE (BE)

dt = 0.05;  % 0.08 limit
h = 0.5;
N = floor((T-t0)/dt);    % # of time intervals
M = floor((b-a)/h);      % # of space intervals
t = linspace(t0,T,N+1)';
x = linspace(a,b,M+1)';

lambda = 0.1;
nmax = 1000;
toll = 1e-5;

% initialize space-time solution matrix
u = zeros(M+1,N+1);

% initial condition
u(:,1) = u0(x);

% temporal loop
for n = 1:N
    % D.B.C.
    u(1,n+1) = 0;
    u(M+1,n+1) = 0;
    % Internal nodes
    [u(2:end-1,n+1),it] = nonlinsolv(u(2:end-1,n),nmax,toll,lambda,dt,h);    % as BE
end

% second lambda
lambda2 = 10;

% initialize space-time solution matrix
u2 = zeros(M+1,N+1);

% initial condition
u2(:,1) = u0(x);

% temporal loop
for n = 1:N
    % D.B.C.
    u2(1,n+1) = 0;
    u2(M+1,n+1) = 0;
    % Internal nodes
    [u2(2:end-1,n+1),it] = nonlinsolv(u2(2:end-1,n),nmax,toll,lambda2,dt,h);    % as BE
end

% mesh
figure(2)
subplot(1, 2, 1);
mesh(t,x,u);
zlim([0 1]);
xlabel('Time');
ylabel('Space of the population');
title('Evolution wave of the propagation of the gene with lambda=0.1');
subplot(1, 2, 2);
mesh(t,x,u2);
zlim([0 1]);
xlabel('Time');
ylabel('Space of the population');
title('Evolution wave of the propagation of the gene with lambda=10');

% lambda is the intensity of selection of the mutant gene or an index of mutation survival advantage
% small lambda --> mutant gene cannot propagate effectively (goes to 0 with time)
% bigger lambda --> mutant gene CAN propagate (almost immediately) --> peak at about t=2.5, then decreases till a regime

%% try with FE

% conditionally stable method: dt <= (h^2)/2

[x_FE, t_FE, u_FE] = FEnonlin(u0,a,b,t0,T,lambda,dt,h);
[x_FE2, t_FE2, u_FE2] = FEnonlin(u0,a,b,t0,T,lambda2,dt,h);

figure(3)
subplot(1, 2, 1);
mesh(t_FE,x_FE,u_FE);
colormap([1 0 0]);
%zlim([0 1]);
xlabel('Time');
ylabel('Space of the population');
%title('Evolution wave of the propagation of the gene with lambda=0.1 and FE');
subplot(1, 2, 2);
mesh(t_FE2,x_FE2,u_FE2);
colormap([1 0 0]);
%zlim([0 1]);
xlabel('Time');
ylabel('Space of the population');
%title('Evolution wave of the propagation of the gene with lambda=10 and FE');

% FEC method for the heat equation: stable if dt <= (h^2)/2 
% in fact: stability if dt<=C*h^2, where C=1/(2*mu), mu=diffusion coeff.=1
% with dt small enough, the same considerations can be done for different values of alpha

%% Propagation of the gene in relation to lambda

% Create a new video file
writerObj = VideoWriter('video_lambda.mov', 'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);

lambda_vect = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10];

figure(4)
for l = lambda_vect
    % initialize space-time solution matrix
    u_l = zeros(M+1,N+1);
    % initial condition
    u_l(:,1) = u0(x);
    % temporal loop
    for n = 1:N
        % D.B.C.
        u_l(1,n+1) = 0;
        u_l(M+1,n+1) = 0;
        % Internal nodes
        [u_l(2:end-1,n+1),it] = nonlinsolv(u_l(2:end-1,n),nmax,toll,l,dt,h);    % as BE
    end
    mesh(t,x,u_l);
    zlim([0 1]);
    title(strcat('lambda = ', num2str(l)));
    frame = getframe(figure(3));
    writeVideo(writerObj,frame);
    pause(0.5);
end

close(writerObj);

% >> lambda, >> 'strenght' of the mutation, >> propagation

%% analysis of the fraction of individuals that have the mutant gene (p)

% integral of p between -2 and 2 in respect to time t

p = zeros(length(t), 1);

for i = 1:length(t)
    p(i) = 0.25 * trapz(u(5:13, i));
    p2(i) = 0.25 * trapz(u2(5:13, i));
end

figure(5)
subplot(1, 2, 1)
plot(t, p, 'b', 'linewidth', 2);
xlabel('Time');
ylabel('p(t)');
title('Evolution of p(t) with lambda=0.1');
subplot(1, 2, 2)
plot(t, p2, 'b', 'linewidth', 2);
xlabel('Time');
ylabel('p(t)');
title('Evolution of p(t) with lambda=10');

sgtitle('Evolution of the fraction of individuals that have the mutant gene - p(t)');

% the two curves start from the same value: p(0) = 0.1875
% with low lambda the gene doesn't propagate and probably tends to 0 (= becomes extinct)
% with high lambda the gene propagates effectively and then reaches a plateau (at around t = 1)

