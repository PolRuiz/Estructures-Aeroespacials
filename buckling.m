function [sig_cr] = buckling(n_d,n_el,x,mat,Tn,Tmat,sig)
%BUCKLING Summary of this function goes here
%   Detailed explanation goes here

x1=zeros(n_el,1);
y1=zeros(n_el,1);
x2=zeros(n_el,1);
y2=zeros(n_el,1);
sig_cr=zeros(n_el,3);


for e=1:n_el
    E=mat(Tmat(e),1);
    A=mat(Tmat(e),2);
    I=mat(Tmat(e),4);
    x1(e) = x(Tn(e,1),1);
    y1(e) = x(Tn(e,1),2);
    z1(e) = x(Tn(e,1),3);
    x2(e) = x(Tn(e,2),1);
    y2(e) = x(Tn(e,2),2);
    z2(e) = x(Tn(e,2),3);
  
    
    L=sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);
    
    sig_cr(e,1)=(pi^2*E*I)/(L^2*A);    
end

sig_cr(:,2)=sig(:,1);

for e=1:n_el
    if sig_cr(e,2)<0
        if sig_cr(e,1)<abs(sig_cr(e,2))
            sig_cr(e,3)=1;
        else
            sig_cr(e,3)=0;
        end        
    end    
    
end

end


