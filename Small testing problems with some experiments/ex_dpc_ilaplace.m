% problem dimension

clear;
n = 60;
tic

% problem ilaplace and its SVD

[A,b,x] = ilaplace(100,1);
[U,S,V] = svd(A);
toc

% different noise levels

b_ex = b;
b_no = randn(100,1);
b_no04 = b_no*1e-04*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
b_no08 = b_no*1e-08*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
b_no14 = b_no*1e-14*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
b04 = b_ex + b_no04;
b08 = b_ex + b_no08;
b14 = b_ex + b_no14;
b04 = b04/sqrt(b04'*b04);
b08 = b08/sqrt(b08'*b08);
b14 = b14/sqrt(b14'*b14);

b = b/sqrt(b'*b);

% plotting results

semilogy(1:100,double(diag(S)),'k-',...
  1:100,abs(U'*b04),'g^',...
  1:100,abs(U'*b08),'bs',...
  1:100,abs(U'*b14),'ro');
  
legend('\sigma_j','|u_j^T b| for \delta_{noise} = 10^{-4}','|u_j^T b| for \delta_{noise} = 10^{-8}','|u_j^T b| for \delta_{noise} = 10^{-14}',4)  
xlabel('singular value number');
%axis([0,400,1e-60,1e0]);  
  
%singvals = double(diag(S2));
%prjctons = abs(double(U2'*b2));
%save('data.mat','singvals','prjctons');