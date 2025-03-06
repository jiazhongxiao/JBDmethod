function relresBound=getrelresBound(target,d)
global B B_bar
%[P,P_bar,~,CA,SL]=gsvd(B,B_bar,0); % In Matlab, gsvd sorts from small to large, svd is the opposite
[P,CA,~]=svd(B(1:d+1,1:d),0); P=P(:,end:-1:1); CA=CA(end:-1:1,end:-1:1);
[P_bar,SL,~]=svd(B_bar(1:d,1:d),0);
if target>0
    i=d-target+1;
    relresBound=abs(B(d+1,d+1)*P(end,i:end)*SL(i:end,i:end)-B_bar(d,d+1)*P_bar(end,i:end)*CA(i:end,i:end));
else
    relresBound=abs(B(d+1,d+1)*P(end,1:-target)*SL(1:-target,1:-target)-B_bar(d,d+1)*P_bar(end,1:-target)*CA(1:-target,1:-target));
end
end