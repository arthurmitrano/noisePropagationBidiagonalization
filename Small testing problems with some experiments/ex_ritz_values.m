% problem dimension

n = 400;

% number of iterations

k = 30;

% problem setting

[A,b,x] = shaw(n);

% SVD of A

[U,Sigma,V] = svd(A);
eigv = (diag(Sigma)).^2;

% Golub-Kahan bidiagonalization

[S,L,W] = bidiag_gk(A,b,k);


for j = [8,18]

  % SVD of L_j

  [P,Theta,Q] = svd(L(1:j,1:j));
  ritz = (diag(Theta)).^2;

  % plotting Ritz values and eigenvalues

  if j == 8,
    subplot(2,1,1);
  else,
    subplot(2,1,2);
  end;
  h = semilogx(ritz,zeros(j,1),'rx',eigv,zeros(n,1),'bo',[1e-40,1e5],[0,0],'k-');
  set(h(1),'LineWidth',2);
  set(h(2),'LineWidth',2);
  set(h(1),'MarkerSize',16);
  set(h(2),'MarkerSize',6);
  set(get(h(1),'Parent'),'FontSize',16);
  if j == 8,
    axis([1e-7,10,-0.5,1]);
  else,
    axis([1e-35,1e5,-0.5,1]);
  end;
%  if j ~= 8,
%    xlabel('Real');
%  end;
%  ylabel('Imag');
  h = legend(sprintf('Ritz values, (\\theta_l^{(%d)})^{ 2}',j),'eigenvalues, \sigma_j^{ 2}');
  set(h,'FontSize',16);
%  print('-dpsc',sprintf('ritz%02d',j));
%  pause;
  
end;