function [X] = Chan_NBS(BSN,BS,R)

%% ��һ��WLS

 %k=X^2+Y^2
for i = 1:BSN                      %BSNΪ��վ����
    k(i) = BS(1,i)^2 + BS(2,i)^2;  %BSΪ��վ����
end

%h = 1/2(Ri^2-ki+k1)
for i =1:BSN-1
    h(i) = 0.5*(R(i)^2 - k(i+1) + k(1));  %ע��k(i+1)
end


%Ga = [Xi,Yi,Ri]
for i = 1:BSN-1
    Ga(i,1) = -BS(1,i+1);
    Ga(i,2) = -BS(2,i+1);
    Ga(i,3) = -R(i);
end

%QΪTDOAϵͳ��Э�������
%Q = cov(R);
Q = eye(BSN-1);

%MS��BS�����Զʱ
za = inv(Ga' * inv(Q) * Ga) * Ga' * inv(Q) * h';

%% �ڶ���WLS
%h'
X1 = BS(1,1);
Y1 = BS(2,1);
h2 = [
    (za(1,1) - X1)^2;
    (za(2,1) - Y1)^2;
     za(3,1)^2
      ];

%Ga'
Ga2 = [1,0;0,1;1,1];

%B'
B2 = [
      za(1,1)-X1,0,0;
      0,za(2,1)-Y1,0;
      0,0,za(3,1)
      ];
  
%za',�����Զʱ
za2 = inv( Ga2' * inv(B2) * Ga' * inv(Q) * Ga * inv(B2) * Ga2) * (Ga2' * inv(B2) * Ga' * inv(Q) * Ga * inv(B2)) * h2;

zp(1,1) = abs(za2(1,1))^0.5 + X1;
zp(2,1) = abs(za2(2,1))^0.5 + Y1;

X = zp;
 

end
