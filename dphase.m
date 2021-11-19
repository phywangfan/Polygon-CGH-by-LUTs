%% ������������˫��λ������Ҷȵ��������һ����һ��Ч��
function [dphs]=dphase(H)
% W ��m*n�ľ���һ��Ϊ1080*1920������Ҫȡ�ڲ���ԭͼ��һ�����ȫϢͼ����



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
