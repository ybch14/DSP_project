function omega=my_esprit(y,n,m)
% y --- signal
% n --- number of sine signals
% m --- the size of convariance matrix
% delta_omega --- the resolution of frequency
N=length(y);
% compute the convariance matrix R
C=xcorr(y,'biased');
R=C(N:N+m-1)';
R=toeplitz(R);
% eigendecomposition, V is eigenvector 
% and D is diag(eigenvalues)
% the eigenvalues from small to big
[V,~]=eig(R);
S=V(:,(m-n+1):(m));
S1=S(1:m-1,:);S2=S(2:m,:);
phi=(S1'*S1)\(S1'*S2);
omega=angle(eig(phi));
return
end