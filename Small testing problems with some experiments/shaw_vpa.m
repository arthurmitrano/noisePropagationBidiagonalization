function [A,b,x] = shaw_vpa(n,digs) 
%SHAW Test problem: one-dimensional image restoration model. 
% 
% [A,b,x] = shaw_vpa(n,digs) 
% 
% Discretization of a first kind Fredholm integral equation with 
% [-pi/2,pi/2] as both integration intervals.  The kernel K and 
% the solution f, which are given by 
%    K(s,t) = (cos(s) + cos(t))*(sin(u)/u)^2 
%    u = pi*(sin(s) + sin(t)) 
%    f(t) = a1*exp(-c1*(t - t1)^2) + a2*exp(-c2*(t - t2)^2) , 
% are discretized by simple quadrature to produce A and x. 
% Then the discrete right-hand b side is produced as b = A*x. 
% 
% The order n must be even. 
%
% This routine uses various precision artihmetic and requires MATLAB 
% Symbolic Toolbox. The parameter "digs" sets the accuracy of Maple's 
% numeric computations, the resultin data [A,b,x] are computed using 
% arthimetic guaranteeing "digs" decimal digits.
%
% See also SHAW, DIGITS, VPA.
 
% Reference: C. B. Shaw, Jr., "Improvements of the resolution of 
% an instrument by numerical solution of an integral equation", 
% J. Math. Anal. Appl. 37 (1972), 83-112. 
 
% Per Christian Hansen, IMM, 08/20/91. 
% Modified by Martin Plesinger, TUL & ASCR, 04/22/09. 
 
% Check input. 
if (rem(n,2)~=0), error('The order n must be even'), end 
 
% Initialization. 
digits(digs);
h = vpa(pi)/n; A = vpa(zeros(n,n)); 
 
% Compute the matrix A. 
co = cos(-vpa(pi)/2 + [.5:n-.5]*h); 
psi = vpa(pi)*sin(-vpa(pi)/2 + [.5:n-.5]*h); 
for i=1:n/2 
  for j=i:n-i 
    ss = psi(i) + psi(j); 
    A(i,j) = ((co(i) + co(j))*sin(ss)/ss)^2; 
    A(n-j+1,n-i+1) = A(i,j); 
  end 
  A(i,n-i+1) = (2*co(i))^2; 
end 
A = A + triu(A,1)'; A = A*h; 
 
% Compute the vectors x and b. 
a1 = 2; c1 = 6; t1 =  .8; 
a2 = 1; c2 = 2; t2 = -.5; 
if (nargout>1) 
  x =   a1*exp(-c1*(-vpa(pi)/2 + [.5:n-.5]'*h - t1).^2) ... 
      + a2*exp(-c2*(-vpa(pi)/2 + [.5:n-.5]'*h - t2).^2); 
  b = A*x; 
end 
