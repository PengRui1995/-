%{
�������ַ����ֱ����BSN=3,4,5,6,7,8,9,10,���㶨λƫ��
����վ���ã�
MS����

%}
clear all
close all 
clc
%BS ����վ
BS=[0,-5,4,-2,7,-7,2,-4,3,1;
    0,8,6,4,3,5,5,2,3,8
]';
%BSN վ��
BSN=10;
%MS ����վ
MS =[8,22];
% expTimes ����ʵ�����
expTimes = 100000;
%MSE �������
MSE =zeros(expTimes,7);
for ind=1:expTimes%������λ
    %R0 ��վ��Դվ�������ֵ
    R0 = zeros(BSN,1);
    for i =1:BSN
        R0(i) = sqrt((BS(i,1)- MS(1,1))^2+(BS(i,2)-MS(1,2))^2);
    end
    %���������õ�����ֵ,
    %R ��������ֵ
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