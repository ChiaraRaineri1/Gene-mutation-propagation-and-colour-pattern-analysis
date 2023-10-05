function [x1,it] = nonlinsolv(x0,nmax,toll,lambda,dt,h)
%% function [x1,it] = nonlinsolv(x0,nmax,toll,lambda,dt,h)
% This routine solves the nonlinear system:
% (I + dt A)*x1 - dt*lambda*x1.*(1-x1) = x0
% where A is the matrix that represents -mu * d^2/d x^2
%
% INPUT: x0 -> right hand side of the non-linear system (used as initial
%              guess as well)
%        nmax -> maximum number of Newton solver method
%        toll -> tollerance treshold for the Newton solver method
%        lambda -> lambda value of the Fisher equation
%        dt -> time-discretization param
%        h -> space discretization param
% OUTPUT: x1 -> solution
%         it -> number of iteration required for convergence

% Initialize iterative method's params
it = 0;
err = toll+1;

% Matrix assembling
M = length(x0);
e = ones(M,1);
B = spdiags([-(dt/h/h)*e, 1 + 2*(dt/h/h)*e, -(dt/h/h)*e],[-1,0,1],M,M); 
% Function definition
f = @(v) B*v - dt*lambda*v.*(1-v) - x0;

% Iterative cycle
while it<nmax && err > toll
   % Compute the jacobian
   J = B - dt*lambda*spdiags([(1-2*x0)],[0],M,M);
   % Newton iteration
   x1 = x0 - J\f(x0);
   % Update of #iteration, error, and solution at previous iteration
   it = it+1;
   err = norm(x1-x0,inf);
   x0 = x1;
end

return