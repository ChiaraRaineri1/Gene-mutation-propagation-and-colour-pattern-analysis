%% FORWARD METHOD FOR NON LINEAR f

function [x, t, u] = FEnonlin(u0,a,b,t0,T,lambda,dt,h)

% Mesh initialization
N = floor((T-t0)/dt);  % # of time intervals  -->  N+1 nodes
M = floor((b-a)/h);    % # of space intervals  -->  M+1 nodes
t = linspace(t0,T,N+1);
x = linspace(a,b,M+1)';

% Initialize space-time solution matrix
u = zeros(M+1,N+1);

% Initial condition
u(:,1) = u0(x);

% Operators definition
e = ones(M-1,1);
A = (spdiags([-e,2*e,-e],[-1,0,1], M-1, M-1))/(h^2);
I = speye(M-1, M-1);

% Temporal loop
for n = 1:N
    % Internal nodes
    u(2:end-1,n+1) = (I - dt*A) * u(2:end-1,n) + lambda * dt * (u(2:end-1,n) - (u(2:end-1,n)).^2);    % usual part + non linear part
    % D.B.C.
    u(1,n+1) = 0;
    u(M+1,n+1) = 0;
end


