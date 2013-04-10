% problem size, number of bidiagonalization iterations, noise level

n = 400;
kmax = 25;
delta_noise = 1e-14;

% problem

[A,b,x] = shaw(n);

% add random noise

b_exact = b;
b_noise = randn(n,1);
b_noise = b_noise*delta_noise*norm(b_exact)/norm(b_noise);
b = b_exact + b_noise;

b_exact = b_exact/norm(b);
b_noise = b_noise/norm(b);
b = b/norm(b);

% bidiagonalization

[S,L,W] = bidiag_gk(A,b,kmax);

% computing the "noise" and "exact" part of the left bidiagonalization vetors

S_noise = [b_noise];
S_exact = [b_exact];

for k = 2:kmax,

  s_noise = - S_noise(:,k-1)*L(k-1,k-1)/L(k,k-1);
  S_noise = [S_noise,s_noise];
  
  s_exact = (A*W(:,k-1) - S_exact(:,k-1)*L(k-1,k-1))/L(k,k-1);
  S_exact = [S_exact,s_exact];
  
end;

% plotting results

subplot(1,2,1);
h = plot(1:n,S_exact(:,18),'k-');
set(h,'LineWidth',2);
axis([0,n,-.2,.2]);
axis square;

subplot(1,2,2);
h = plot(1:n,S_noise(:,18),'k-');
set(h,'LineWidth',2);
axis([0,n,-.2,.2]);
axis square;
