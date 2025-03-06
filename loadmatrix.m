function [m,n,p]=loadmatrix(matrixpath,csvpath)
% load matrices A and L and save the informations to csvpath (if provided)
global A L
S=load(matrixpath);
A=S.Problem.A;
if size(A,1)>size(A,2)
    A=A';
end
m=size(A,1); n=size(A,2);
if isfield(S.Problem,'L')
    L=S.Problem.L;
    n=size(L,1);
else
    i=[1:n,2:n,1:n-1];
    j=[1:n,1:n-1,2:n];
    v=[3*ones(1,n),ones(1,2*n-2)];
    p=n;
    L=sparse(i,j,v,p,n);
end
if csvpath~=""
    matrixname=erase(matrixpath,"../IRJBDmatrices/");
    matrixname=erase(matrixname,".mat");
    kappa=svds(A,1)/svds(A,1,'smallest');
    variablenames={'A','m','n','p','nnz','kappa'};
    T=table(matrixname,m,n,p,nnz(A)+nnz(L),kappa,'VariableNames',variablenames);
    writetable(T,csvpath,'WriteMode','append');
end
end