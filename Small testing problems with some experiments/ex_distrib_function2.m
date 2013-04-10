
% problem setting

n = 400;
ks = [1:10];
delta_no = 1e-14;

%[A,b,x] = shaw(n);

%A = rand(n);
%b = rand(n,1);

b = rand(n,1);
A = gallery('fiedler',b);   

% adding the noise

b_ex = b;
b_no = randn(n,1);
b_no = b_no*delta_no*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
b = b_ex + b_no;
b_norm = sqrt(b'*b);
%x = x/b_norm;
b = b/b_norm;
b_ex = b_ex/b_norm;
b_no = b_no/b_norm;

% for different values of k 

for k = ks,

  % SVD of A and bidiagonalization of A

  [U,Sigma,V] = svd(A);
  [S,L,W]     = bidiag_gk(A,b,k,0,0,-1,2);

  % Riemann-Stieltjes distribution function omega 

  nodes = diag(Sigma);
  nodes = nodes.^2; 
  nodes = nodes(n:-1:1);

  nodes2 = [];
  for j = 1:n,
    nodes2 = [nodes2; nodes(j); nodes(j)];
  end;
  nodes2 = [nodes2; 1e10];

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
    weight_sums2 = [weight_sums2; weight_sums(j); weight_sums(j)];
  end;

  % Ritz values and the approximation of the distribution function

  ritz_nodes = [];
  ritz_weight_sums = [];
  [P,Theta,Q] = svd(L);
  ritz_nodes = diag(Theta).^2;
  ritz_weights = P(1,:).^2;

  ritz_weight_sums = [];
  for j = 1:k,
    ritz_weight_sums(j) = sum(ritz_weights(j:k));
  end;
  
  ritz_nodes2 = [];
  ritz_nodes2 = [ritz_nodes2; 1e10];
  for j = 1:k,
    ritz_nodes2 = [ritz_nodes2; ritz_nodes(j); ritz_nodes(j)];
  end;

  ritz_weight_sums2 = [];
  for j = 1:k,
    ritz_weight_sums2 = [ritz_weight_sums2; ritz_weight_sums(j); ritz_weight_sums(j)];
  end;
  ritz_weight_sums2 = [ritz_weight_sums2; 1e-100];
   
  % plotting graphs

  %H = loglog(ritz_nodes2,ritz_weight_sums2,'r-',nodes2,weight_sums2,'k-',ritz_nodes2,ritz_weight_sums2,'rx');
  H = loglog(ritz_nodes2,ritz_weight_sums2,'r-',nodes2,weight_sums2,'k-');
  set(get(H(1),'Parent'),'FontSize',16);
  set(H(2),'LineWidth',2);
  set(H(1),'LineWidth',3);
  %set(H(3),'MarkerSize',16);

  xlabel('nodes');
  ylabel('cumulated weights');
  G = legend('approximed distribution function','original distribution function',2);
  G = legend('approximed distribution function \omega^{(k)}','original distribution function \omega',2);
  
  % suitable for shaw(400) problem
  % axis([1e-16,1e2,1e-16,1e1]);
 
  % suitable for rand(400) problem
  % axis([1e-1,1e5,1e-2,1e1]);
 
  % suitable for fiedler problem
  axis([1e-2,1e6,1e-5,1e1]);
  
  %print('-dpsc','interlacing.eps');

  pause(1);

end;


