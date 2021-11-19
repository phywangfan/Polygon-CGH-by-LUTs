%% 待做工作：将双相位和振幅灰度调整结合在一起试一试效果
function [dphs]=dphase(H)
% W 是m*n的矩阵，一般为1080*1920，现在要取内部与原图像一样大的全息图出来



[m,n]=size(H);
ag=angle(H);  % phase of hologram
A=abs(H)/max(max(abs(H))); % amplitude of hologram
pa=ag+acos(A/2);  
pa(pa<0)=pa(pa<0)+2*pi;
pb=ag-acos(A/2);
pb(pb<0)=pb(pb<0)+2*pi;

dphs=zeros(m,n);
for k=1:m  % row
    if(mod(k,2)==1)
        dphs(k,1:2:n)=pa(k,1:2:n);
        dphs(k,2:2:n)=pb(k,2:2:n); 
    else
        dphs(k,2:2:n)=pa(k,2:2:n);
        dphs(k,1:2:n)=pb(k,1:2:n);
    end
end


dphs=(dphs-min(dphs(:)))./(max(dphs(:))-min(dphs(:)));
end
% imwrite(flipud(phs),'phs_BeforeRangeAdjust.bmp')
