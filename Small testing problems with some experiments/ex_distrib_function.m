% It generates approximation for the distribution function for several problems

for PROBLEM = [1:5],

  % problem dimensions, number of bidiagonalization steps, noise levels
  
  if PROBLEM == 1,
    n = 400;
    k = 25;
    delta_nos = [1e-14,1e-10,1e-6,1e-4,1e-2];
  elseif PROBLEM == 2,
    n = 100;
    k = 30;
    delta_nos = [1e-13,1e-10,1e-7,1e-2,1e-1];
  elseif PROBLEM == 3,
    n = 400;
    k = 30;
    delta_nos = [1e-14,1e-13,1e-10,1e-6,1e-7,1e-4,1e-2,1e-1];
  elseif PROBLEM == 4,
    n = 400;
    k = 200;
    delta_nos = [1e-14,1e-13,1e-10,1e-6,1e-7,1e-4,1e-2,1e-1];
  elseif PROBLEM == 5,
    n = 400;
    k = 100;
    delta_nos = [1e-14,1e-13,1e-10,1e-6,1e-7,1e-4,1e-2,1e-1];
  end; 

  for delta_no = delta_nos,

    % problem setting

    if      PROBLEM == 1, [A,b,x] = shaw(n);
    elseif  PROBLEM == 2, [A,b,x] = ilaplace(n,1);
    elseif  PROBLEM == 3, [A,b,x] = phillips(n);
    elseif  PROBLEM == 4, [A,b,x] = heat(n,1);
    elseif  PROBLEM == 5, [A,b,x] = heat(n,3);
    end;
    b_ex = b;
    b_no = randn(n,1);
    b_no = b_no*delta_no*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
    b = b_ex + b_no;
    b_norm = sqrt(b'*b);
    x = x/b_norm;
    b = b/b_norm;
    b_ex = b_ex/b_norm;
    b_no = b_no/b_norm;

    % SVD of A and bidiagonalization of A

    [U,Sigma,V] = svd(A);
    [S,L,W]     = bidiag_gk(A,b,k,0,0,-1,2);

    % Riemann-Stieltjes distribution function omega 

    nodes = diag(Sigma);
    nodes = nodes.^2; 
    nodes = nodes(n:-1:1);

    nodes2 = [];
    for j = 1:n,
      nodes2 = [nodes2;nodes(j);nodes(j)];
    end;
    nodes2 = [nodes2;1e10];

    weights = U'*b;
    weights = weights.^2;
    weights = weights(n:-1:1);

    weight_sums = [];
    for j = 1:n,
      weight_sums(j) = sum(weights(1:j));
    end;

    weight_sums2 = [];
    weight_sums2 = [weight_sums2; 1e-100];
    for j = 1:n,
      weight_sums2 = [weight_sums2;weight_sums(j);weight_sums(j)];
    end;

    % Ritz values (note that the singular values and Ritz values have different ordering)

    ritz_nodes = [];
    ritz_weights = [];
    for j = 1:k,
      [P,Theta,Q] = svd(L(1:j,1:j));
      ritz_nodes = [ritz_nodes;Theta(j,j)^2];
      ritz_weights = [ritz_weights;P(1,j)^2];
    end;

    % plotting graphs

    H = loglog(nodes2,weight_sums2,'k-',ritz_nodes,ritz_weights,'ro',[1e-100,1e+100],[1,1]*delta_no^2,'k-.');
    set(get(H(1),'Parent'),'FontSize',16);
    set(H(1),'LineWidth',2);
    set(H(2),'LineWidth',2);
    set(H(2),'MarkerSize',6);
    set(H(2),'MarkerFaceColor','red');

    xlabel('nodes');
    ylabel('cumulated weights');
    G = legend('\omega(\lambda): nodes \sigma_j^2, weights |(b/\beta_1,u_j)|^2',...
      'nodes (\theta_1^{(k)})^2 and its weights |(p_1^{(k)},e_1)|^2',...
      '\delta_{noise}^2',2);
  
    if      PROBLEM == 1, axis([1e-41,1e2,delta_no^2.214285,1e1]);
    elseif  PROBLEM == 2, axis([1e-41,1e2,delta_no^2.214285,1e1]);
    elseif  PROBLEM == 3, axis([1e-21,1e2,0.001*delta_no^2.214285,1e1]);
    elseif  PROBLEM == 4, axis([1e-81,1e2,0.001*delta_no^2.214285,1e1]);
    elseif  PROBLEM == 5, axis([1e-41,1e2,0.001*delta_no^2.214285,1e1]);
    end;

    % printing figures

    if      PROBLEM == 1, print('-dpsc',sprintf('omegafce_shaw(%d)_%e.eps',n,delta_no));
    elseif  PROBLEM == 2, print('-dpsc',sprintf('omegafce_ilaplace(%d,1)_%e.eps',n,delta_no));
    elseif  PROBLEM == 3, print('-dpsc',sprintf('omegafce_phillips(%d)_%e.eps',n,delta_no));
    elseif  PROBLEM == 4, print('-dpsc',sprintf('omegafce_heat(%d,1)_%e.eps',n,delta_no));
    elseif  PROBLEM == 5, print('-dpsc',sprintf('omegafce_heat(%d,3)_%e.eps',n,delta_no));
    end;

  end;
end;

