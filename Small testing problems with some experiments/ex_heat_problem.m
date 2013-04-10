% noise level, number of iterations

noise_level = 1e-3;
k = 100;

% problem size

n = 400;

% values of parameter kappa in the heat problem

for kappa = 1:0.2:5,

  % problem 
  
  [A,b,x] = heat(n,kappa);

  % random noise

  b_exact = b;
  b_noise = randn(n,1);
  b_noise = b_noise*noise_level*norm(b_exact)/norm(b_noise);

  b = b_exact + b_noise;
  b = b/norm(b);

  % bidiagonalization with and without reorthogonalization

  [Sy,Ly,Wy] = bidiag_gk(A,b,k);
  [Sn,Ln,Wn] = bidiag_gk(A,b,k,0,0,0,0);

  % weight corresponding to the smallest Ritz value with (y) and without (n)
  % reorthogonalization
  
  datay = [];
  datan = [];
  for j = 1:k,
    [Uy,Sigmay,Vy] = svd(Ly(1:j,1:j));
    [Un,Sigman,Vn] = svd(Ln(1:j,1:j));
    datay = [datay, Uy(1,j)];
    datan = [datan, Un(1,j)];
  end;

  % SVD of A

  [UU,SSigma,VV] = svd(A);

  % plotting the discrete Picard condition

  subplot(1,2,1)
  h = semilogy(1:n,diag(SSigma),'r-',1:n,abs(UU'*b),'b.');
  legend('\sigma_j','|(u_j,b)|');
  axis([0,n,1e-10,1])

  % plotting the smalest Ritz values

  subplot(1,2,2)
  h = semilogy([1,k],noise_level*[1,1],'k--',1:k,abs(datay),'r-',1:k,abs(datan),'b-.');
  %set(h(1),'LineWidth',2);
  %set(h(2),'LineWidth',2);
  %set(h(2),'MarkerSize',12);
  %set(h(3),'LineWidth',2);
  %set(h(3),'MarkerSize',12);
  legend('noise level',...
    sprintf('heat(400,%1.1f), reort. YES',kappa),...
    sprintf('heat(400,%1.1f), reort. NO',kappa));
  axis([0,k,1e-4,1]);

  print('-dpsc',sprintf('heat_kappa=%1.1f.eps',kappa));

end;
