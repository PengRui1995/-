
clear all 
close all
clc
%{
	��վ��λ����������
%}

BSN =3;
BS = [0,0;20,20;0,20;20,0];
BS = 50*BS(1:BSN,:);
MS = [50,79];
R = zeros(BSN-1,1);
for i = 1:BSN-1
	R(i)=sqrt((MS(1)-BS(i+1,1))^2+(MS(2)-BS(i+1,2))^2)-sqrt((MS(1)-BS(1,1))^2+(MS(2)-BS(1,2))^2)+0*randn(1);
end
%% ��λ���򲿷�
Pos = Chan_3BS(BS,R)%δ�򻯰汾
Pos = Chan1(BS,R)% �򻯰汾
%%
%Pos = Chan1(BS,R);