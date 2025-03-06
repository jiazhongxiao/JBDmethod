function [d,relresBoundVec,flag]=IRJBD(target,k,adjust,reorth,tol,maxit,lsqrtol,lsqrmaxit)
d=0; l=abs(target);
relresBoundVec=zeros(maxit,1);
for iter=1:maxit
    [d,relresBound,flag]=JBD(target,d,k,reorth,tol,lsqrtol,lsqrmaxit);
    relresBoundVec(iter)=max(relresBound);
    if flag==0
        break
    end
    fprintf("iter=%d, relresBound=%e\n",iter,relresBoundVec(iter));
    [lambda2,mu2]=getShifts(sign(target)*(l+adjust),d,relresBound);
    imRestartJBD(l+adjust,k,lambda2,mu2);
    d=l+adjust;
    if iter==1 || mod(iter,10)==0
        testJBD(d);
    end
end
relresBoundVec=relresBoundVec(1:iter);
end