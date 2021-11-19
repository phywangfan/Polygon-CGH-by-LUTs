function [F] = ifc(F,u,v)
ee=1e-10;
p=logical(abs(v)<ee);q=logical(abs(u)<ee);r=logical(abs(v+u)<ee);
if any(p(:)) || any(q(:)) || any(r(:))  %% 表明p q r 只要任意一个非空就可以算下去
    if  any(p(:)) && any(abs(u(p))<ee)       %  1) 当u=0 且v=0时
        k=logical(p==1&q==1); % 找到v=0且u=0的索引
        F(k)=0.5;
    elseif any(abs(v(q))>ee)                  %  2) 当u=0 v!=0时
        k=logical(p~=1&q==1); % 找到v!=0且u=0的索
        F(k)=(1-exp(-2*pi*1j*v(k)))./((2*pi*v(k)).^2)-1j./(2*pi*v(k));
    elseif any(abs(u(p))>ee)                  %  3) 当 u!=0 v=0时
        k=logical(q~=1 & p==1); % 找到v=0且u!=0的索
        F(k)=(exp(-2*pi*1j*u(k))-1)./((2*pi*u(k)).^2)+1j*exp(-2*pi*1j*u(k))./(2*pi*u(k));
    elseif any(r(:)) && any(abs(v(r))>ee)     %  4) 当u+v=0,但u!=0，v!=0
        k=logical(r==1&p~=1); % 找到v!=0且u!=0的索
        F(k)=(1-exp(2*pi*1j*v(k)))./((2*pi*v(k)).^2)+1j./(2*pi*v(k));
    end
end

