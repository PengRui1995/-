function Pos =Chan2(BSN,BS,R)
%函数的功能：基于Chan氏算法N站定位的函数
%函数的描述：输入多个台站的横纵坐标，已经测量到个各站到信号源站的距离差，就可以估计出信号源站的位置
%函数的使用：Pos=Chan2(BSN,BS,R)
%输入：
%     BSN:用于定位的台站数
%     BS:三个台站的横纵坐标,每行代表一个接收站坐标
%     R:测量到个各站到信号源站的距离差
%输出：
%     Pos:估计出信号源站的位置的横纵坐标,会存在四组坐标，需要利用先验知识进行取舍。
%例子：无;
%注意事项：利用函数的适用范围,在已知信号源离接收站较远的条件下使用,但实验计算发现近站条件下亦可用。
%文档日期：~
%标签：~
%创建日期：~
%最后更新日期：~
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
%Q 测量误差的协方差矩阵
%Q = cov(R);
Q = eye(BSN-1);
%% MS与BS较远时
%z 第一次WPS对z的估计值
z= inv(G'*inv(Q)*G)*G'*inv(Q)*h;
%第二次WPS
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
Pos =zp';

%% 消除模糊定位 
%计算后者中离第一次估计点最近的位置
z_t= ones(4,2);
z_t(:,1)=z_t(:,1)*z(1);
z_t(:,2)=z_t(:,2)*z(2);
d= Pos-z_t;
for i =1:4
dis(i) = d(i,1)^2+d(i,2)^2;
end
[val,ind]=min(dis);
Pos = Pos(ind,:);
end