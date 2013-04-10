clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting of two problems 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem shaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sh_A,sh_b,sh_x] = shaw(400);

% random noise with different levels

sh_n = randn(400,1);
sh_n = sh_n/norm(sh_n);

sh_n1 = sh_n*1e-12; 
d1 = norm(sh_n1)/norm(sh_b); % delta 2.14e-14
sh_n2 = sh_n*1e-8;  
d2 = norm(sh_n2)/norm(sh_b); % delta 2.14e-10
sh_n3 = sh_n*1e-4;
d3 = norm(sh_n3)/norm(sh_b); % delta 2.14e-6
sh_n4 = sh_n*1e-2;
d4 = norm(sh_n4)/norm(sh_b); % delta 2.14e-4
sh_n5 = sh_n*1e-0;
d5 = norm(sh_n5)/norm(sh_b); % delta 2.14e-2

sh_noise = [sh_n1,sh_n2,sh_n3,sh_n4,sh_n5];
sh_d = [d1,d2,d3,d4,d5];

% number bidiagonalization iterations

sh_maxit = [20,18,13,10,6];

% vertical axis range

sh_ymax = [1e+02,1e+02,1e+01,1e+01,1e+01];
sh_ymin = [1e-14,1e-14,1e-07,1e-06,1e-03];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem ilaplace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[il_A,il_b,il_x,il_t] = ilaplace(100,1);

% random noise with different levels

il_n = randn(100,1);
il_n = il_n/norm(il_n);

il_n1 = il_n*1e-12; 
d1 = norm(il_n1)/norm(il_b); % delta 2.14e-13
il_n2 = il_n*1e-9; 
d2 = norm(il_n2)/norm(il_b); % delta 2.14e-10
il_n3 = il_n*1e-6; 
d3 = norm(il_n3)/norm(il_b); % delta 2.14e-7
il_n4 = il_n*1e-1; 
d4 = norm(il_n4)/norm(il_b); % delta 2.14e-2
il_n5 = il_n*1e-0; 
d5 = norm(il_n5)/norm(il_b); % delta 2.14e-1

il_noise = [il_n1,il_n2,il_n3,il_n4,il_n5];
il_d = [d1,d2,d3,d4,d5];

% number bidiagonalization iterations

il_maxit = [25,23,19,9,5];

% vertical axis range

il_ymax = [1e+02,1e+02,1e+01,1e+01,1e+01];
il_ymin = [1e-14,1e-10,1e-08,1e-03,1e-02];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bidiagonalization, plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem shaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:5,

  % without reorthogonalization

  [u0,l0,v0] = bidiag_gk(sh_A,sh_b+sh_noise(:,i),sh_maxit(i),1,eps,0,0);
  alfa = diag(l0);
  beta = diag(l0,-1);
  ratio = beta./alfa;
  for j=2:sh_maxit(i),
    ratio(j) = ratio(j)*ratio(j-1);
  end;    
  
  % plotting results
  
%  subplot(1,2,1);
  h = semilogy(1:sh_maxit(i),alfa,'-.',1:sh_maxit(i),beta,'--',1:sh_maxit(i),ratio,'-');
  set(h,'LineWidth',2);
  set(get(h(1),'Parent'),'YMinorTick','off');
  axis([1,sh_maxit(i),sh_ymin(i),sh_ymax(i)]);
%  legend('\alpha_1,...,\alpha_{k_{noise}+2}','\beta_2,...,\beta_{k_{noise}+3}','estimate (3.7)',3);
  legend('\alpha_{j-1}','\beta_j','estimate \rho_k',3);
%  title(sprintf('Bidiagonal elements of L, SHAW(400), \\delta_{noise} = %e without reorthogonalization',sh_d(i)));
  print('-dpsc',sprintf('./shaw(400)_delta=%e_without.eps',sh_d(i)));

  % with reorthogonalization

  [u2,l2,v2] = bidiag_gk(sh_A,sh_b+sh_noise(:,i),sh_maxit(i),1,eps,-1,2);
  alfa = diag(l2);
  beta = diag(l2,-1);
  ratio = beta./alfa;
  for j=2:sh_maxit(i),
    ratio(j) = ratio(j)*ratio(j-1);
  end;    

  % plotting results

%  subplot(1,2,2);
  h = semilogy(1:sh_maxit(i),alfa,'-.',1:sh_maxit(i),beta,'--',1:sh_maxit(i),ratio,'-');
  set(h,'LineWidth',2);
  set(get(h(1),'Parent'),'YMinorTick','off');
  axis([1,sh_maxit(i),sh_ymin(i),sh_ymax(i)]);
%  legend('\alpha_1,...,\alpha_{k_{noise}+2}','\beta_2,...,\beta_{k_{noise}+3}','estimate (3.7)',3);
  legend('\alpha_{j-1}','\beta_j','estimate \rho_k',3);
%  title(sprintf('Bidiagonal elements of L, SHAW(400), \\delta_{noise} = %e with reorthogonalization',sh_d(i)));
  print('-dpsc',sprintf('./shaw(400)_delta=%e_with.eps',sh_d(i)));
%  print('-dpsc',sprintf('./shaw(400)_delta=%e.eps',sh_d(i)));

end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem ilaplace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:5,

  % without reorthogonalization
  
  [u0,l0,v0] = bidiag_gk(il_A,il_b+il_noise(:,i),il_maxit(i),1,eps,0,0);
  alfa = diag(l0);
  beta = diag(l0,-1);
  ratio = beta./alfa;
  for j=2:il_maxit(i),
    ratio(j) = ratio(j)*ratio(j-1);
  end;    

  % plotting results

%  subplot(1,2,1);
  h = semilogy(1:il_maxit(i),alfa,'-.',1:il_maxit(i),beta,'--',1:il_maxit(i),ratio,'-');
  set(h,'LineWidth',2);
  set(get(h(1),'Parent'),'YMinorTick','off');
  axis([1,il_maxit(i),il_ymin(i),il_ymax(i)]);
%  legend('\alpha_1,...,\alpha_{k_{noise}+2}','\beta_2,...,\beta_{k_{noise}+3}','estimate (3.7)',3);
  legend('\alpha_{j-1}','\beta_j','estimate \rho_k',3);
%  title(sprintf('Bidiagonal elements of L, ILAPLACE(100,1), \\delta_{noise} = %e without reorthogonalization',il_d(i)));
  print('-dpsc',sprintf('./ilaplace(100,1)_delta=%e_without.eps',il_d(i)));

  % with reorthogonalization

  [u2,l2,v2] = bidiag_gk(il_A,il_b+il_noise(:,i),il_maxit(i),1,eps,-1,2);
  alfa = diag(l2);
  beta = diag(l2,-1);
  ratio = beta./alfa;
  for j=2:il_maxit(i),
    ratio(j) = ratio(j)*ratio(j-1);
  end;    

  % plotting results

%  subplot(1,2,2);
  h = semilogy(1:il_maxit(i),alfa,'-.',1:il_maxit(i),beta,'--',1:il_maxit(i),ratio,'-');
  set(h,'LineWidth',2);
  set(get(h(1),'Parent'),'YMinorTick','off');
  axis([1,il_maxit(i),il_ymin(i),il_ymax(i)]);
%  legend('\alpha_1,...,\alpha_{k_{noise}+2}','\beta_2,...,\beta_{k_{noise}+3}','estimate (3.7)',3);
  legend('\alpha_{j-1}','\beta_j','estimate \rho_k',3);
%  title(sprintf('Bidiagonal elements of L, ILAPLACE(100,1), \\delta_{noise} = %e with reorthogonalization',il_d(i)));
  print('-dpsc',sprintf('./ilaplace(100,1)_delta=%e_with.eps',il_d(i)));
%  print('-dpsc',sprintf('./ilaplace(100,1)_delta=%e.eps',il_d(i)));

end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

