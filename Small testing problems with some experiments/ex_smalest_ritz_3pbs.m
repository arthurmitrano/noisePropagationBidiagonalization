% comparison of three problems

clear;

% plots ragnes

axsizes = [...
  0,100,1e-18,1;...
  0,100,1e-16,1;...
  0,100,1e-14,1;...
  0,100,1e-12,1;...
  0,50,1e-10,1;...
  0,50,1e-8,1;...
  0,50,1e-6,1;...
  0,50,1e-5,1;...
  0,50,1e-4,1;...
  0,50,1e-2,1];

% noise levels 

levels = [-12:2:-8,-7:-1];

for jj = 1:10,

  % choose the noise level

  noise_level = 10^levels(jj)
  
  % with(out) reorthogonalization

  for reort = [0,1],

    k = 100;

    % problem shaw
  
    n = 400;
    [A,b,x] = shaw(n);
    b_exact = b;
    b_noise = randn(n,1);
    b_noise = b_noise*noise_level*norm(b_exact)/norm(b_noise);
    b = b_exact + b_noise;

    % bidiagonalization with and without reorthogonalization

    if reort == 0,
      [S,L,W] = bidiag_gk(A,b,k);
    else,
      [S,L,W] = bidiag_gk(A,b,k,0,0,0,0);
    end;

    % weigths corresponding to the smalest Ritz values

    data1 = [];
    for j = 1:k,
      [U,Sigma,V] = svd(L(1:j,1:j));
      data1 = [data1, U(1,j)];
    end;

    % problem ilaplace

    n = 100;
    [A,b,x] = ilaplace(n,1);
    b_exact = b;
    b_noise = randn(n,1);   
    b_noise = b_noise*noise_level*norm(b_exact)/norm(b_noise);  
    b = b_exact + b_noise;

    % bidiagonalization with and without reorthogonalization

    if reort == 0,
      [S,L,W] = bidiag_gk(A,b,k);
    else,
      [S,L,W] = bidiag_gk(A,b,k,0,0,0,0);
    end;

    % weigths corresponding to the smalest Ritz values

    data2 = [];
    for j = 1:k,
      [U,Sigma,V] = svd(L(1:j,1:j));
      data2 = [data2, U(1,j)];
    end;

    % problem phillips

    n = 400;
    [A,b,x] = phillips(n);
    b_exact = b;
    b_noise = randn(n,1);
    b_noise = b_noise*noise_level*norm(b_exact)/norm(b_noise);
    b = b_exact + b_noise;

    % bidiagonalization with and without reorthogonalization

    if reort == 0,
      [S,L,W] = bidiag_gk(A,b,k);
    else,
      [S,L,W] = bidiag_gk(A,b,k,0,0,0,0);
    end;

    % weigths corresponding to the smalest Ritz values

    data3 = [];
    for j = 1:k,
      [U,Sigma,V] = svd(L(1:j,1:j));
      data3 = [data3, U(1,j)];
    end;

    % plotting results

    h = semilogy([1,k],noise_level*[1,1],'k--',1:k,abs(data1),'r-x',1:k,abs(data2),'g-o',1:k,abs(data3),'b-^');
    %set(h(1),'LineWidth',2);
    %set(h(2),'LineWidth',2);
    %set(h(2),'MarkerSize',12);
    %set(h(3),'LineWidth',2);
    %set(h(3),'MarkerSize',12);
    %set(h(4),'LineWidth',2);
    %set(h(4),'MarkerSize',12);
    legend('noise level','shaw(400)','ilaplace(100,1)','phillips(400)')

    axis(axsizes(jj,:));
    if reort == 0,
      title(sprintf('with reorthogonalization, noise level \\delta = 10^{%d}',levels(jj)));
      print('-dpsc',sprintf('threepbs_reorty_noiselev-%02d',abs(levels(jj))));
    else,
      title(sprintf('without reorthogonalization, noise level \\delta = 10^{%d}',levels(jj)));
      print('-dpsc',sprintf('threepbs_reortn_noiselev-%02d',abs(levels(jj))));
    end;

  end;
end;

