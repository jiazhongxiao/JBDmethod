function [lambda2,mu2]=getShifts(l,d,relresBound)
% abs(l): dimension after restart，d: current dimension
global B B_bar
%[P,P_bar,~,CA,SL]=gsvd(B,B_bar,0); % In Matlab, gsvd sorts from small to large, svd is the opposite
c=svd(B(1:d+1,1:d),0);
c=c(end:-1:1);
s=svd(B_bar(1:d,1:d),0);
if l>0
    lambda2=c(1:d-l).^2;
    relgap=abs(c(d-l+1)-relresBound(1)-c(1:d-l))/c(d-l+1);
    for i=1:d-l
        if relgap(i)<1e-3
            continue;
            lambda2(i)=0.0;
        end
    end
    mu2=1-lambda2;
else
    lambda2=c(d:-1:1-l).^2;
    relgap=abs(c(-l)-relresBound(-l-3)-c(d:-1:1-l))/c(-l);
    for i=1:d+l
        if relgap(i)<1e-3
            lambda2(i)=1.0;
        end
    end
    mu2=1-lambda2;
end
end