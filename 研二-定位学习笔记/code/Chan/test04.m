%{
利用三种方法分别计算BSN=3,4,5,6,7,8,9,10,计算定位偏差
接收站设置：
MS设置

%}
clear all
close all 
clc
%BS 接收站
BS=[0,-5,4,-2,7,-7,2,-4,3,1;
    0,8,6,4,3,5,5,2,3,8
]';
%BSN 站数
BSN=10;
%MS 发射站
MS =[8,22];
% expTimes 独立实验次数
expTimes = 100000;
%MSE 估计误差
MSE =zeros(expTimes,7);
for ind=1:expTimes%独立定位
    %R0 各站到源站距离的真值
    R0 = zeros(BSN,1);
    for i =1:BSN
        R0(i) = sqrt((BS(i,1)- MS(1,1))^2+(BS(i,2)-MS(1,2))^2);
    end
    %加入噪声得到测量值,
    %R 距离差测量值
    R = zeros(BSN-1,1);
    noise =0.01;
    for i=1:BSN-1
        R(i) = R0(i+1)-R0(1)+noise*randn(1);
    end
    for i=4:10
        Pos = Chan3(i,BS,R(1:i-1));
        MSE(ind,i-3) = (Pos(1)-MS(1))^2+(Pos(2)-MS(2))^2;
    end
    
end
MSE_mean =mean(MSE);