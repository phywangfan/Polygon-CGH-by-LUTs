% --------- Running this procedure at first to generate look-up tables,
% --------- then running 'Main1....m' procedure to generate CGH. 
% --------- By Fan Wang, Chiba University, Japan, 2021/11/19.
% --------- wangfan@chiba-u.jp, phy_wangfan@163.com
%% fundamental parameters
clear
Lx=15.36; % size of hologram; also size of object in x and y axis
Ly=8.64;
depth=20;   % depth in z-axis direction
lamda=532e-6;
k=2*pi/lamda;
smp_x=1920;   % number of samples
smp_y=1080;
dp_x=Lx/(smp_x-1);
dp_y=Ly/(smp_y-1);          % the period of pixel
d=500;%max(Lx,Ly)*dp_x/lamda*2; 
Lzo=d+depth;

mag_xy=round(log10(max(Lx,Ly)),0);   % order of magnitudes of x-,and y-axis direction
mag_z=round(log10(depth),0);             % order of magnitudes of z-axis direction                              
prcs_xy=10^mag_xy/1e3;       % the precision of x,y responding to its order of magnitudes
prcs_z=0.04;%depth/1e3*2;
%% define the original frequency coordinates
prcs_uv=1e-2;
dp_u0=round(1/Lx,abs(log10(prcs_uv)));   % precision of frequency u0
dp_v0=round(1/Ly,abs(log10(prcs_uv))); 
u0=-dp_u0/2-(smp_x/2-1)*dp_u0:dp_u0:dp_u0/2+(smp_x/2-1)*dp_u0;
v0=(dp_v0/2+(smp_y/2-1)*dp_v0:-dp_v0:-dp_v0/2-(smp_y/2-1)*dp_v0)';
%% e^(j*2π*(u'+v')) and  e^(j*2π*u') 
fx1=(floor(-Lx/10/prcs_xy)*prcs_xy:prcs_xy:ceil(Lx/10/prcs_xy)*prcs_xy)';  % factor of Horizontal, (a11+a12)=(x2-x1) or a11=x2-x1
fy1=(floor(-Ly/10/prcs_xy)*prcs_xy:prcs_xy:ceil(Ly/10/prcs_xy)*prcs_xy)';   % a21+a22  or a21
fz1=(floor(-depth/10/prcs_z)*prcs_z:prcs_z:ceil(depth/10/prcs_z)*prcs_z);  % a31+a32 or a13 =z2-z1
fz1(fz1==0)=1;
LU1=zeros(length(fx1)*length(fz1),smp_x);
LV1=zeros(smp_y,length(fy1)*length(fz1));
for i=1:length(fz1)
    LU1(length(fx1)*(i-1)+1:length(fx1)*i,:)=exp(1j*2*pi*lamda/2*(u0-fx1./fz1(i)/lamda).^2.*fz1(i));   
    LV1(:,length(fy1)*(i-1)+1:length(fy1)*i)=exp(1j*2*pi*lamda/2*(v0-fy1'./fz1(i)/lamda).^2.*fz1(i));
end
len_fx1=length(fx1);
len_fy1=length(fy1);
%% e^(j*2π*(a14*u0+a24*v0+a34*w0)) 
a14=[-Lx/2,Lx/2];  % a14∈[-Lx/2,Lx/2] [-7.5, 7.5]
a24=[-Ly/2,Ly/2];  % a24∈[-Ly/2,Ly/2] [-4, 4] 
a34=[d,(d+depth)]; % a34∈[d+depth/2,d-depth/2]  [414, 464]
A=[min(a14)/max(a34),max(a14)/min(a34)]/lamda;             % min_A=min(a14)/max(a34)
A=max(abs(A));
B=[min(a24)/max(a34),max(a24)/min(a34)]/lamda;
B=max(abs(B));
temp1=u0(1)-prcs_uv:-prcs_uv:u0(1)-ceil(A/prcs_uv)*prcs_uv;     % possible u0 value in the negative direction
temp2=u0(end)+prcs_uv:prcs_uv:u0(end)+ceil(A/prcs_uv)*prcs_uv;  % possible u0 value in the positive direction
ua=[sort(temp1),u0(1):prcs_uv:u0(end),temp2];   % possible value of (u0-a14/a34/lamda)
clear temp1 temp2
temp1=v0(1)+prcs_uv:prcs_uv:v0(1)+ceil(B/prcs_uv)*prcs_uv;      % possible v0 value in the positive direction
temp2=v0(end)-prcs_uv:-prcs_uv:v0(end)-ceil(B/prcs_uv)*prcs_uv; % possible v0 value in the negative direction
vb=[sort(temp1,'descend'),v0(1):-prcs_uv:v0(end),temp2]';  % possible value of (v0-a24/a34/lamda)
a34=ceil(a34(1)/prcs_z)*prcs_z:prcs_z:ceil(a34(2)/prcs_z)*prcs_z;  % possible value of a34
LU0=zeros(length(a34),length(ua));
LV0=zeros(length(vb),length(a34));
for i=1:length(a34)
    LU0(i,:)=exp(1j*2*pi*lamda/2*ua.^2*a34(i));  
    LV0(:,i)=exp(1j*2*pi*lamda/2*vb.^2*a34(i));  
end
SMP_y=length(v0(1):-prcs_uv:v0(end));
SMP_x=length(u0(1):prcs_uv:u0(end));
%% e^(-j2π*(a11+a12)*u0) && e^(-j2π*(a21+a22)*v0)
LU2=zeros(length(fx1),smp_x);
LV2=zeros(smp_y,length(fy1));
for i=1:length(fx1)
    LU2(i,:)=exp(-1j*2*pi*fx1(i)*u0); 
end
for i=1:length(fy1)
    LV2(:,i)=exp(-1j*2*pi*fy1(i)*v0);
end
%% PCA compressed for LU0 and LV0
meanLU0=mean(LU0);
[coeff,score,~,~,explained] = pca(LU0);
explained(:,2)=cumsum(explained);
avi=find((100-explained(:,2))<1e-13, 1 );
coeff_LU0=coeff(:,1:avi);
score_LU0=score(:,1:avi);
clear  coeff score explained avi

meanLV0=mean(LV0');
[coeff,score,~,~,explained] = pca(LV0');
explained(:,2)=cumsum(explained);
avi=find((100-explained(:,2))<1e-13, 1 );
coeff_LV0=coeff(:,1:avi);
score_LV0=score(:,1:avi);
clear  coeff score explained avi
%% PCA compressed for LU1 and LV1
meanLU1=mean(LU1);
[coeff,score,~,~,explained] = pca(LU1);
explained(:,2)=cumsum(explained);
avi=find((100-explained(:,2))<1e-13, 1 );
coeff_LU1=coeff(:,1:avi);
score_LU1=score(:,1:avi);
clear  coeff score explained avi

meanLV1=mean(LV1');
[coeff,score,~,~,explained] = pca(LV1');
explained(:,2)=cumsum(explained);
avi=find((100-explained(:,2))<1e-13, 1 );
coeff_LV1=coeff(:,1:avi);
score_LV1=score(:,1:avi);
clear  coeff score explained avi

%%
parameter.lamda=lamda;
parameter.smp_x=smp_x;
parameter.smp_y=smp_y;
parameter.Lx=Lx; parameter.Ly=Ly;
parameter.delta=prcs_xy;
parameter.d=d;
parameter.prcs_xy=prcs_xy;
parameter.prcs_z=prcs_z;
parameter.depth=depth;
parameter.dp_x=dp_x;
parameter.dp_y=dp_y;
parameter.prcs_uv=prcs_uv;
parameter.SMP_y=SMP_y;
parameter.SMP_x=SMP_x;
parameter.len_fx1=len_fx1;
parameter.len_fy1=len_fy1;
parameter.eps=eps;

%% save data to hardware
fname='E:\hologramLUTs\parameter_plg.mat';
save(fname,'parameter');
fname='E:\hologramLUTs\coeff_LU0.mat';
save(fname,'coeff_LU0');
fname='E:\hologramLUTs\coeff_LU1.mat';
save(fname,'coeff_LU1');
fname='E:\hologramLUTs\coeff_LV0.mat';
save(fname,'coeff_LV0');
fname='E:\hologramLUTs\coeff_LV1.mat';
save(fname,'coeff_LV1');

fname='E:\hologramLUTs\score_LU0.mat';
save(fname,'score_LU0');
fname='E:\hologramLUTs\score_LU1.mat';
save(fname,'score_LU1');
fname='E:\hologramLUTs\score_LV0.mat';
save(fname,'score_LV0');
fname='E:\hologramLUTs\score_LV1.mat';
save(fname,'score_LV1');
fname='E:\hologramLUTs\meanLU0.mat';
save(fname,'meanLU0');
fname='E:\hologramLUTs\meanLU1.mat';
save(fname,'meanLU1');
fname='E:\hologramLUTs\meanLV0.mat';
save(fname,'meanLV0');
fname='E:\hologramLUTs\meanLV1.mat';
save(fname,'meanLV1');

fname='E:\hologramLUTs\LU2.mat';
save(fname,'LU2');
fname='E:\hologramLUTs\LV2.mat';
save(fname,'LV2');



