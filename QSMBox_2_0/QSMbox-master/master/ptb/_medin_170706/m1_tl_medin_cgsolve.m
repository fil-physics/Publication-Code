function [x, rsd, iter] = m1_tl_medin_cgsolve(A, b, tol, maxiter, verbose, x0)
% m1_tl_medin_cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = m1_tl_medin_cgsolve(A, b, tol, maxiter, verbose, x0)
%
% A - Either an NxN matrix, or a function handle.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    If and integer greater than 0, print out progress every 'verbose' iters.
%
% x0 - initial solution
%
% Adapted from original code by Justin Romberg
% Email: jrom@acm.caltech.edu
% Created: October 2005
%
% Modified for QSM by Ildar Khalidov, WCMC
% Modified by Tian Liu on 2011.02.02
% Last modified by Julio Acosta-Cabronero

disp(' '); disp('CG solver')

c=0; % JAC
matrix_size=size(b);
b=b(:);

if (nargin < 6)
    x = zeros(length(b),1);
else 
    x=x0(:);
end

if (nargin < 5), verbose = 1; end

implicit = isa(A,'function_handle');

if (nargin < 6)
    r = b;
else
    if (implicit), r = b - A(reshape(x,matrix_size)); r=r(:);  
    else r = b - A*x;  end
end

d = r;
delta = r'*r;
delta0 = b'*b;
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0); 

while numiter<maxiter && delta>tol^2*delta0 
  if implicit; q = A(reshape(d,matrix_size)); q=q(:);  
  else q = A*d;  end
 
  alpha = delta/(d'*q);
  x = x + alpha*d;
  
  if (mod(numiter+1,50) == 0)
    if (implicit), r = b - reshape(A(reshape(x,matrix_size)),size(b));  
    else r = b - A*x;  end
  else
    r = r - alpha*q;
  end
  
  deltaold = delta;
  delta = r'*r;
  beta = delta/deltaold;
  d = r + beta*d;
  numiter = numiter + 1; c=c+1; % JAC
  
  if sqrt(delta/delta0)<bestres
    bestx = x;
    bestres = sqrt(delta/delta0);
  end
  
    verbose=2; % JAC
    printprogress = false;
    if verbose==1 && mod(numiter,verbose)==0
        printprogress = true;
    elseif verbose==2 && mod(numiter/5,verbose)==0
        printprogress = true;
    end
    if printprogress
        disp(sprintf('#%d. Residual = %.3e (tol=%.0e, max_iter=%d)',numiter,bestres,tol,maxiter));
    end
end

verbose=1; % JAC
if verbose
    disp(' ')
    disp(sprintf('CG stop @%d. Final residual = %.3e',numiter,sqrt(delta/delta0)));
end

x = reshape(bestx,matrix_size);
rsd = bestres;
iter = numiter;
