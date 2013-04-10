clear;

% problem formulation

problem_size = 400;
delta_no = 1e-14;
max_it = 25;

reort = -1;
refin = 2;

[A,b,x] = shaw(problem_size);
b_ex = b;
b_no = randn(problem_size,1);
b_no = b_no*delta_no*sqrt(b_ex'*b_ex)/sqrt(b_no'*b_no);
b = b_ex + b_no;
b_norm = sqrt(b'*b);
x = x/b_norm;
b = b/b_norm;
b_ex = b_ex/b_norm;
b_no = b_no/b_norm;

% bidiagonalization and SVD of the bidiagonal matrix

[S,L,W] = bidiag_gk(A,b,max_it,0,0,reort,refin);
[U,sigma,V] = svd(A);

% computing sizes of components of left bidiagonalization vectors S_{j-1}, S_j
% in the signal space span([U_1,...,U_j]) and 
% in the noise space span([U_{j+1},...,U_n])

for k = 1:max_it,
  if k == 1,
    sspace = [ 0 , norm( U(:,1:k)' * S(:,k) ) ];
    nspace = [ 0 , norm( U(:,k+1:problem_size)' * S(:,k) ) ];
  else,
    sspace = [ sspace ; norm( U(:,1:k)' * S(:,k-1) ) , norm( U(:,1:k)' * S(:,k) ) ];
    nspace = [ nspace ; norm( U(:,k+1:problem_size)' * S(:,k-1) ) , norm( U(:,k+1:problem_size)' * S(:,k) ) ];
  end;
end;

% plotting figures

for k = 1:max_it,
  if k > 1,
    subplot(5,5,k) 
    h = plot(sin(0:0.01:pi/2),cos(0:0.01:pi/2),'k-',...
      sspace(k,1),nspace(k,1),'b^',sspace(k,2),nspace(k,2),'ro');
    title(sprintf('s_{%d} -> s_{%d}',k-1,k))
    axis square
    axis([0,1,0,1])
    set(get(h(1),'Parent'),'XTick',[]);
    set(get(h(1),'Parent'),'YTick',[]);
    set(h(2),'MarkerFaceColor',[0 0 1]);
    set(h(3),'MarkerFaceColor',[1 0 0]);
  end;
end;





return;

% print -dpsc figure1.eps
pause

data = [2,6,7,11,14,16];
for k = 1:6,
    subplot(1,6,k) 
    h = plot(sin(0:0.01:pi/2),cos(0:0.01:pi/2),'k-',...
      sspace(data(k),1),nspace(data(k),1),'b^',sspace(data(k),2),nspace(data(k),2),'ro');
    title(sprintf('s_{%d} -> s_{%d}',data(k)-1,data(k)))
    axis square
    axis([0,1,0,1])
    set(get(h(1),'Parent'),'XTick',[]);
    set(get(h(1),'Parent'),'YTick',[]);
    set(h(2),'MarkerFaceColor',[0 0 1]);
    set(h(3),'MarkerFaceColor',[1 0 0]);
end;

% print -dpsc figure2.eps
pause;

data = [17,18,19,20,21,22];
for k = 1:6,
    subplot(1,6,k) 
    h = plot(sin(0:0.01:pi/2),cos(0:0.01:pi/2),'k-',...
      sspace(data(k),1),nspace(data(k),1),'b^',sspace(data(k),2),nspace(data(k),2),'ro');
    title(sprintf('s_{%d} -> s_{%d}',data(k)-1,data(k)))
    axis square
    axis([0,1,0,1])
    set(get(h(1),'Parent'),'XTick',[]);
    set(get(h(1),'Parent'),'YTick',[]);
    set(h(2),'MarkerFaceColor',[0 0 1]);
    set(h(3),'MarkerFaceColor',[1 0 0]);
end;

% print -dpsc figure3.eps