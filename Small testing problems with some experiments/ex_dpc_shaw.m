% problem dimension

clear;
n = 60;
digs = 110;

% approximation of problem shaw (computed using standard arithmetic)

[A,b,x] = shaw(n);
[U,S,V] = svd(A);

% problem shaw (coputed using various precision arithmetic)

[Avpa,bvpa,xvpa] = shaw_vpa(n,digs);
[Uvpa,Svpa,Vvpa] = svd(Avpa);

% plotting results

h = semilogy(1:n,diag(S),'k-',1:n,double(diag(Svpa)),'r--',1:n,abs(double(Uvpa'*bvpa)),'bo');
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
set(h(3),'LineWidth',2);
%legend('\sigma_j(A), [A,b,x] = shaw(60)','\sigma_j(A), [A,b,x] = shaw\_vpa(60)','| U^T b^{exact} | of A, [A,b,x] = shaw\_vpa(60)',3);
legend('\sigma_j(A), shaw(60)','\sigma_j(A), shaw\_vpa(60)','| U^T b^{exact} |, shaw\_vpa(60)',3);
print -dpsc shaw_vpa_DPC.eps