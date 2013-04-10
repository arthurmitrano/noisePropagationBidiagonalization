% problem and its singular value decomposition

[A,b,x] = shaw(400);
[U,Sigma,V] = svd(A);

% plotting and printing selected left singular vectors 

for j = 1:5,

  subplot(1,5,j), plot(1:400,U(:,j),'k');
  axis([0,400,-0.1,0.1]);
  axis square;
  title(sprintf('u_{%d}',j));
  
end;

print -dpsc uvect1.eps
pause

for j = 6:10,

  subplot(1,5,j-5), plot(1:400,U(:,j),'k');
  axis([0,400,-0.1,0.1]);
  axis square;
  title(sprintf('u_{%d}',j));
  
end;

print -dpsc uvect2.eps
pause
