function c = CT(x1,x2,x3,lamda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c1=exp(-1j*2*pi/lamda/2*(x1^2+x2^2)/x3);
c2=exp(-1j*2*pi/lamda*x3);
c=c1*c2;
end

