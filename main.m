addpath('IRJBDtools');
matrixlist=dir("../IRJBDmatrices/*.mat");
matrixnames=[];
for i=1:size(matrixlist,1)
    name=string(matrixlist(i).name);
    name=erase(name,".mat");
    matrixnames=[matrixnames,name];
end
global A L U U_hat V_prime B B_bar;
write=1; % Whether to write results into a csv file
targets=[5,-5];
ks=[25,50];
adjust=3;
tol=1e-8;
maxit=1000;
reorth=1; % 0: no reorthogonalization，1: full reorthogonalization 2: reorthogonalization only on V_prime and U, not on U_hat
lsqrtol=10*eps;
for matrixname=matrixnames
    [m,n,p]=loadmatrix(sprintf("../IRJBDmatrices/%s.mat",matrixname),"");
    lsqrmaxit=10*n;
    global A L U U_hat V_prime B B_bar;
    Rbound=sqrt(norm(A,1)*norm(A,inf)+norm(L,1)*norm(L,inf));
    rng(2024); % random seed
    u1=normalize(randn(m,1),"norm");
    for target=targets
        for k=ks
            fprintf("matrixname: %s, target: %d, k: %d\n",matrixname,target,k);
            U=zeros(m,k+1); U(:,1)=u1; U_hat=zeros(p,k); V_prime=zeros(m+p,k+1);
            B=zeros(k+1,k+1); B_bar=zeros(k,k+1);
            tic;
            [d,relresBoundVec,flag]=IRJBD(target,k,adjust,reorth,tol,maxit,lsqrtol,lsqrmaxit);
            time=toc;

            [CA,SL,Y,Z,X]=computeGSVD(target,d,lsqrtol,lsqrmaxit);
            relres=0;
            for i=1:abs(target)
                r1=norm(A*X(:,i)-CA(i,i)*Y(:,i));
                r2=norm(L*X(:,i)-SL(i,i)*Z(:,i));
                r3=norm(SL(i,i)*A'*Y(:,i)-CA(i,i)*L'*Z(:,i));
                relres=max(relres,sqrt(r1^2+r2^2+r3^2)/Rbound);
            end

            if write==1
                variablenames={'A','target','k','IRJBDiter','TRJBDiter','IRJBDtime','Res_b','IRJBDRes','TRJBDRes','invB','invB_bar','flag'};
                T=table(matrixname,target,k,size(relresBoundVec,1),0,time,relresBoundVec(end),relres,0,...,
                    norm(inv(B(1:d,1:d))),norm(inv(B_bar(1:d,1:d))),flag,'VariableNames',variablenames);
                writetable(T,'../results.csv','WriteMode','append');
            end
            %save(sprintf("../relresBoundVec/%s,%d,%d,%d.mat",matrixname,target,k,size(relresBoundVec,1)),"relresBoundVec");
        end
    end
end
rmpath('IRJBDtools');