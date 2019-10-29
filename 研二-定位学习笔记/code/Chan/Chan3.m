function  Pos =Chan3(BSN,BS,R)
%函数的功能：基于Chan氏算法N站定位的函数,考虑近站条件下对源位置做进一步的估计
%函数的描述：输入多个台站的横纵坐标，已经测量到个各站到信号源站的距离差，就可以估计出信号源站的位置
%函数的使用：Pos=Chan3(BSN,BS,R)
%输入：
%     BSN:用于定位的台站数
%     BS:三个台站的横纵坐标,每行代表一个接收站坐标
%     R:测量到个各站到信号源站的距离差
%输出：
%     Pos:估计出信号源站的位置的横纵坐标,会存在四组坐标，需要利用先验知识进行取舍。
%例子：无;
%注意事项：利用函数的适用范围,在已知信号源离接收站较远的条件下使用,但实验计算发现近站条件下亦可用。
%符号说明:公式中出现了h，h'，G,G'，z，z’表示不同矩阵名称，这里改用ha,hb,Ga,Gb,za,zb作为标识符，以免与其转置矩阵混淆。

%k  = xi^2+yi^2
k =zeros(BSN,1);
for i =1:BSN
    k(i) = BS(i,1)^2+BS(i,2)^2;
end
% ha =0.5*(ri1^2+k1) i=2,..m
ha = zeros(BSN-1,1);
for i =1:BSN-1
    ha(i) = 0.5*(R(i)^2-k(i+1)+k(1));
end
%Ga = -[xi1,yi1,ri1]
Ga = zeros(BSN-1,3);
for i=1:BSN-1
    Ga(i,1:2)=BS(i+1,:)-BS(1,:);
    Ga(i,3) = R(i);
end
Ga =-Ga;
%Q 测量误差的协方差矩阵
%Q =cov(R);
Q = 0.5*(eye(BSN-1)+ones(BSN-1,BSN-1));
%% 第一次WLS远站情况下
za = (Ga'/Q*Ga)\Ga'/Q*ha;

%% 第一次WLS近站情况下
%Ba =diag（r2,r3..rm）ri=ri1+r1;
r =zeros(BSN-1,1);
r = R+za(3)*ones(BSN-1,1);
Ba =diag(r');
%
Sigma = Ba*Q*Ba;
W = inv(Sigma);
% za 重新估计
za =(Ga'*W*Ga)\Ga'*W*ha;

%% 第二次WLS近站情况下
%hb 
hb =zeros(3,1);
hb(1) =(za(1)-BS(1,1))^2;
hb(2) =(za(2)-BS(1,2))^2;
hb(3) =za(3)^2;
%Gb 
Gb=[1,0;0,1;1,1];
%Bb
Bb =diag(za-[BS(1,:),0]');
%zCov =(Ga'*W*Ga)^-1
zCov =inv(Ga'*W*Ga);
%SigmaB  =4*Bb*zCov*Bb
SigmaB = 4*Bb*zCov*Bb;
%zb
zb=(Gb'/SigmaB*Gb)\Gb'/SigmaB*hb;
%%
res = sqrt(zb);
%% Pos
Pos(1,:) = [res(1) res(2)]+BS(1,:);
Pos(2,:) = [res(1) -res(2)]+BS(1,:);
Pos(3,:) = [-res(1) res(2)]+BS(1,:);
Pos(4,:) = [-res(1) -res(2)]+BS(1,:);

%%  选择与第一次wls的结果相近的点
dis= zeros(4,1);
for i=1:4
    dis(i) = [Pos(i,:)-za(1:2,1)']* [Pos(i,:)-za(1:2,1)']';
end
[val,ind] = min(dis);
%%
Pos = Pos(ind,:);

end