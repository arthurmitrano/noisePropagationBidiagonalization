% problem dimensions, number of iterations, noise level

n = 400;
k = 25;
delta_no = 1e-14;

% problem setting, adding the noise

[A,b,x]     = shaw(n);
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

%SJ   = [];
%SJP1 = [];

% plotting the computed bidiagonalization vectors

for j = 1:k-1,
%for j = [7,12,16],

  % bidiagonalization vectors in the SVD bases, \alpha_j, and \beta_{j+1}

  sj   = U'*S(:,j);
  sjp1 = U'*S(:,j+1);
  wj   = V'*W(:,j);
  aj   = L(j,j);
  bjp1 = L(j+1,j);

  % print iteration step and the \beta{j+1}/\alpha_j ratio

  j, bjp1/aj

  %SJ   = [SJ,   log(abs(aj*sj))    ];
  %SJP1 = [SJP1, log(abs(bjp1*sjp1))];

  % plotting the computed vectors

  h = semilogy(1:n,abs(Sigma*wj),'k:o',1:n,abs(aj*sj),'b--+',1:n,abs(bjp1*sjp1),'r-');
  set(get(h(1),'Parent'),'FontSize',16);
  set(h(1),'LineWidth',2);
  %set(h(1),'MarkerSize',12);
  set(h(2),'LineWidth',2);
  set(h(2),'MarkerSize',12);
  set(h(3),'LineWidth',2);
  %set(h(3),'MarkerSize',12);
  h = get(1,'Children');
  set(h(1),'FontSize',16);
  legend(sprintf('| \\Sigma V^T w_{%d} |',j),...
    sprintf('| \\alpha_{%d} U^T s_{%d} |',j,j),...
    sprintf('| \\beta_{%d} U^T s_{%d} |',j+1,j+1));
  
  %if j == 7,
    axis([1,30,1e-30,1]);
  %else,
  %  axis([1,25,1e-25,1e-5]);
  %end;
  %title(sprintf('\\beta_{j+1} s_{j+1} = A w_j - \\alpha_j s_j, j = %d, vse v bazi levych singularnich vektoru A',j));

  % printing figures

  print('-dpsc',sprintf('rovnice_%02d.eps',j));
  pause;

end;

%surf(1:29,[1:40,40:-1:1],[SJ(1:40,:);SJP1(40:-1:1,:)],[SJ(1:40,:)-SJP1(1:40,:);SJP1(40:-1:1,:)-SJ(40:-1:1,:)]);














