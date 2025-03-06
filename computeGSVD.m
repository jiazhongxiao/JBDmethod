function [CA,SL,Y,Z,X]=computeGSVD(l,d,lsqrtol,lsqrmaxit)
global A L U U_hat V_prime B B_bar;
[P,P_bar,W,CA,SL]=gsvd(B(1:d+1,1:d),B_bar(1:d,1:d),0);
X=zeros(size(A,2),abs(l));
if l>0
    CA=CA(d:-1:d-l+1,d:-1:d-l+1);
    SL=SL(d:-1:d-l+1,d:-1:d-l+1);
    Y=U(:,1:d+1)*P(:,d:-1:d-l+1);
    Z=U_hat(:,1:d)*P_bar(:,d:-1:d-l+1);
    for i=1:l
        [X(:,i),~]=lsqr([A;L],V_prime(:,1:d)*W(:,d-i+1),lsqrtol,lsqrmaxit);
    end
else
    CA=CA(1:-l,1:-l);
    SL=SL(1:-l,1:-l);
    Y=U(:,1:d+1)*P(:,1:-l);
    Z=U_hat(:,1:d)*P_bar(:,1:-l);
    for i=1:-l
        [X(:,i),~]=lsqr([A;L],V_prime(:,1:d)*W(:,i),lsqrtol,lsqrmaxit);
    end
end
end