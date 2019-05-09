function [X,D]=jeig(A,B);
L=chol(B,'lower');
G=inv(L);
C=G*A*G';
[Q,D]=schur(C);
X=G'*Q;