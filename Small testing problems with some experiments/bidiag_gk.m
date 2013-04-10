function [U,L,V] = bidiag_gk(A,s,k,p,tol,reort,refin)

% [U,L,V] = bidiag_gk(A,s) or L = bidiag_gk(A,s)   ... mandatory inputs, it 
%           computes either all three factors or only the bidiagonal matrix
% [U,L,V] = bidiag_gk(A,s,k)                       ... number of iterations
% [U,L,V] = bidiag_gk(A,s,k,p)                     ... square/rectangular L
% [U,L,V] = bidiag_gk(A,s,k,p,tol)                 ... warning threshold 
% [U,L,V] = bidiag_gk(A,s,k,p,tol,reort)           ... reorthogonalization
% [U,L,V] = bidiag_gk(A,s,k,p,tol,reort,refin)     ... iterative refinement
% =========================================================================
% Input arguments:
% A     ... matrix,
% s     ... starting nonzero vector (e.g. the right-hand side in Ax=b),
% k     ... number of iterations of the bidiagonalization process
%           (each iteration consists of the left and right half-step),
% p     ... proces stops with square or rectangular bidiagonal matrix:
%           if p == 0 or p == false (or p not present), then
%
%    A^T * U_k = V_k * L_k^T,    where U_k = U(:,1:k), V_k = V(:,1:k) and
%
%          L_k = / a_1          0  \      
%                | b_2 a_2         |          
%                |     ... ...     |
%                \  0      b_k a_k /
%
%    is square lower bidiagonal matrix, if p == 1 or p == true or p > 1, then
%
%      A * V_k = U_{k+1} * L_{k+},    where
%
%       L_{k+} = / a_1          0      \ = / L_k             \
%                | b_2 a_2             |   \ b_{k+1} * e_k^T /
%                |     ... ...         |
%                |         b_k a_k     |
%                \  0          b_{k+1} /
%
%    is rectangular lower bidiagonal matrix,
%
% tol   ... threshold/tolerance, default value is 2^-52 = 2.2204*10^-16
%           (if a_j < tol or b_j < tol for some j, then algorithm warns),  
%           THIS FUNCTION IS TEMPORARILY DISABLED!!!
% reort ... reorthogonalization
%           ==  0         without reorthogononalization,
%           ==  r > 0     reorthogonalization against last r vectors,
%           == -1         full reorthogonalization (against all vecs. default),
% refin ... iterative refinement (multiple reorthogonalization),
%           ==  0         without reortogonalization,
%           ==  1         only one reorthogonalization, CGS-like,
%           ==  2         double reorthogonalization, ICGS-like (default),
%           ==  t         t-times reorth
% (reort and refin default values are recommended.)
%
% =========================================================================
%
% [Golub, Kahan: Calculating the singular values and pseudo-inverse of 
% a matrix, SIAM J. Numer. Anal., Ser. B, 2 (1965), pp. 205--224]

% Martin Plesinger, 2005--2009, TU Liberec & AS CR


% checking inputs 

if nargin < 2, 
    L = [];
    U = [];
    V = [];
    fprintf('Invalid input parameters.\n');
    return; 
elseif nargin == 2,
    k = max(size(A));
    p = 0; 
    tol = 2^-52;
    reort = -1;
    refin = 2;
elseif nargin == 3, 
    p = 0; 
    tol = 2^-52;
    reort = -1;
    refin = 2;
elseif nargin == 4,
    tol = 2^-52;
    reort = -1;
    refin = 2;
elseif nargin == 5,
    reort = -1;
    refin = 2;
elseif nargin == 6,
    refin = 2;
end;

% checking outputs

if (nargout == 2) | (nargout > 3),
    L = [];
    U = [];
    V = [];
    fprintf('Invalid output parameters, can be only L or [U,L,V].\n');
    return; 
end;    

% problem size

[n,m] = size(A);

% correction of the number of iterations (no more its. than the problem size)

if m < n,
    if k > m,   k = m;      p = true;   end;
elseif m == n
    if k > m,   k = m;      end;
    if k == m,  p = false;  end;
elseif m > n,
    if k > n,   k = n;      end;
    if k == n,  p = false;  end;
end;

% print information about the bidiagonalization

fprintf('\n');
fprintf('bidiagonalization technique: Golub-Kahan (Lanczos-like) technique\n');
fprintf('matrix dimensions   = [%d, %d]\n',size(A,1),size(A,2));
fprintf('iterations count    = %d',k);
if p, 
    fprintf('.5'); 
end; 
fprintf('\n');
fprintf('tolerance           = %e\n',tol);
fprintf('reorthogonalization = ');
if reort < 0,
    fprintf('full');
elseif reort == 0,
    fprintf('off');
else
    fprintf('over last %d vectors',reort);
end;
fprintf('\n');
fprintf('iteration refining  = %dx\n',refin);

tic

% nomalization of starting vector, i.e. computing of vector u_1

u1 = s;
beta = sqrt(u1'*u1);
u1 = u1/beta;

% vector v_1

v1 = A'*u1;
alfa = sqrt(v1'*v1);
v1 = v1/alfa;

% update matrices L, U, V

L(1,1) = alfa;
U = u1;
V = v1;
%if alfa < tol, fprintf('Info: alfa_1 = %e is smaller than tolerance.\n',alfa); end;

% the main cycle

for i = 2:k,

    % vector u_i

    u2 = A*v1 - u1*alfa;
    for l = 1:refin,
        
        % multiple (default double) reortogonalization, ICSG-like iterative refinement
        
        for j = i-1:-1:((reort<0)*1 + (reort>=0)*max(1,i-reort)),
            
            % full reorthogonalization
            % if reort=-1 ... 1 + 0*max = 1
            %    (i-1):(-1):1
            %
            % no reorthogonalization
            % if reort= 0 ... 0 + 1*max(1,i-0) = i
            %    (i-1):(-1):i
            %
            % reorthogonalization over last k vectors
            % if reort= k ... 0 + 1*max(1,i-k) = 1 [if i<=k], i-k [if i>k]
            %    (i-1):(-1):1     [if i<=k]
            %    (i-1):(-1):(i-k) [if i>k]
            
            u2 = u2 - (U(:,j)'*u2)*U(:,j);
        end;
    end;
    beta = sqrt(u2'*u2);
    u2 = u2/beta;

    % update matrices L, U
    
    L(i,i-1) = beta;
    U = [U,u2];
    %if beta < tol, fprintf('Info: beta_%d = %e is smaller than tolerance.\n',i,beta); end;

    % vector v_i
    
    v2 = A'*u2 - v1*beta;
    for l = 1:refin,
        for j = i-1:-1:((reort<0)*1 + (reort>=0)*max(1,i-reort)),
          v2 = v2 - (V(:,j)'*v2)*V(:,j);
        end;
    end;
    alfa = sqrt(v2'*v2);
    v2 = v2/alfa;

    % update matrices L, V
    
    L(i,i) = alfa;
    V = [V,v2];
    %if alfa < tol, fprintf('Info: alfa_%d = %e is smaller than tolerance.\n',i,alfa); end;

    % rename vectors
    
    u1 = u2;
    v1 = v2;

end;

% computing the last half-step if a rectangular bidiagonal matrix is wanted

if p,

    % vector u_{k+1}

    u2 = A*v1 - u1*alfa;
    for j =  k:-1:(reort<0)*1 + (reort>=0)*max(1,k+1-reort)
      u2 = u2 - (U(:,j)'*u2)*U(:,j);
    end;
    beta = sqrt(u2'*u2);
    u2 = u2/beta;

    % update matrices L, U

    L(k+1,k) = beta;
    U = [U,u2];
    %if beta < tol, fprintf('Info: beta_%d = %e is smaller than tolerance.\n',k+1,beta); end;

end;    

% print the CPU-time

toc

% checking outputs

if nargout <= 1,
    U = L;
end;

return;
