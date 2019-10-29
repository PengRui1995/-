%{
��BSN>3ʱ�Ķ�λ�㷨���в��ԣ�����Ч������д�ɺ���Chan2.
ʱ�䣺2019-10-8
%}
clear all
close all
clc
BSN =4;
BS = [0,0;sqrt(3),0;0.5*sqrt(3),1.5;-0.5*sqrt(3),1.5];
BS = 50*BS(1:BSN,:);
MS = [1000,1000];
%R0 ��վ��Դվ�������ֵ
R0 = zeros(BSN,1);
for i =1:BSN
    R0(i) = sqrt((BS(i,1)- MS(1,1))^2+(BS(i,2)-MS(1,2))^2);
end
%���������õ�����ֵ,
%R ��������ֵ
R = zeros(BSN-1,1);
noise =0.001;
for i=1:BSN-1
    R(i) = R0(i+1)-R0(1)+noise*randn(1);
    
end
%% ��λ����
PosS = Chan_NBS(BSN,BS',R)
PosL = Chan2(BSN,BS,R)
PosM = Chan3(BSN,BS,R)
%k ki= xi^2+yi^2;
K =zeros(BSN,1);
for i=1: BSN
    K(i) = BS(i,1)^2+BS(i,2)^2;
end
%h =0.5*(ri2-ki+k1)
h=zeros(BSN-1,1);
for i=1:BSN-1
    h(i)=0.5*(R(i)^2-K(i+1)+K(1));
end
%G = -[Xi1 yi1 ri1] 
G = zeros(BSN-1,3);
for i=1:BSN-1
    G(i,1:2) = BS(i+1,:)-BS(1,:);
    G(i,3) = R(i);
end
G = -G;
%Q ��������Э�������
Q = cov(R);%������
Q = eye(BSN-1);

%% MS��BS��Զʱ
%z ��һ��WPS��z�Ĺ���ֵ
z= inv(G'*inv(Q)*G)*G'*inv(Q)*h;
%�ڶ���WPS
%h2 
h2= zeros(3,1);
h2(1)=(z(1)-BS(1,1))^2;
h2(2)=(z(2)-BS(1,2))^2;
h2(3)=z(3)^2;
%G2 
G2=[1,0;0,1;1,1];
B2=diag([z(1)-BS(1,1),z(2)-BS(1,2),z(3)]);

%% z2 =[x,y] 
z2 = inv(G2'*inv(B2)*G'*inv(Q)*G*inv(B2)*G2)*G2'*inv(B2)*G'*inv(Q)*G*inv(B2)*h2;
res = sqrt(abs(z2));
zp(:,1) = [res(1);res(2)]+BS(1,:)';
zp(:,2) = [res(1);-res(2)]+BS(1,:)';
zp(:,3) = [-res(1);res(2)]+BS(1,:)';
zp(:,4) = [-res(1);-res(2)]+BS(1,:)';
Pos =zp'
