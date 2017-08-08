function omega=my_root_music(y,n,m)
% y is the input signal
% n is the number of single-frequency signal
% m is the size of covariance matrix
N=length(y);
% compute the convariance matrix R
C=xcorr(y,'biased');
R=C(N:N+m-1)';
R=toeplitz(R);
% eigendecomposition, V is eigenvector 
% and D is diag(eigenvalues)
% the eigenvalues from small to big
[V,~]=eig(R);
G=V(:,1:(m-n));
C=G*G';
a=zeros(1,m);
for i=1:m
   a(m-i+1)=sum(diag(C,i-1)); 
end
ra=roots(a);
[~,I]=sort(abs(abs(ra)-1));
omega=angle(ra(I(1:n)));
end