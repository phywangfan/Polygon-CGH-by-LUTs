function [F] = ifc(F,u,v)
ee=1e-10;
p=logical(abs(v)<ee);q=logical(abs(u)<ee);r=logical(abs(v+u)<ee);
if any(p(:)) || any(q(:)) || any(r(:))  
    if  any(p(:)) && any(abs(u(p))<ee)       
        k=logical(p==1&q==1);  
        F(k)=0.5;
    elseif any(abs(v(q))>ee)                 
        k=logical(p~=1&q==1); 
        F(k)=(1-exp(-2*pi*1j*v(k)))./((2*pi*v(k)).^2)-1j./(2*pi*v(k));
    elseif any(abs(u(p))>ee)                
        k=logical(q~=1 & p==1); 
        F(k)=(exp(-2*pi*1j*u(k))-1)./((2*pi*u(k)).^2)+1j*exp(-2*pi*1j*u(k))./(2*pi*u(k));
    elseif any(r(:)) && any(abs(v(r))>ee)
        k=logical(r==1&p~=1); 
        F(k)=(1-exp(2*pi*1j*v(k)))./((2*pi*v(k)).^2)+1j./(2*pi*v(k));
    end
end

