%Fang定位算法
clc;
clear;
C = 300000000; % 光速(米每秒)
R=2000;       %小区半径
% 蜂窝系统各基站坐标
X1 = 0;
Y1 = 0;
X2 = R*(sqrt(3));
Y2 = 0;
X3 = R*(sqrt(3)); 
Y3 = R*(-3/2);
X4 = R*(-sqrt(3)/2);
Y4 = R*(-3/2);
X5 = R*(-sqrt(3));
Y5 = 0;
X6 = R*(-sqrt(3)/2);
Y6 = R*(3/2);
X7 = R*(sqrt(3)/2);
Y7 = R*(3/2);

% KM = XM^2 + YM^2
K1 = Ka(X1,Y1);
K2 = Ka(X2,Y2);
K3 = Ka(X3,Y3);
K4 = Ka(X4,Y4);
K5 = Ka(X5,Y5);
K6 = Ka(X6,Y6);
K7 = Ka(X7,Y7);

% 随机产生MS的位置（x，y）
 u = rand(1);
 y = (R/2)*(1-sqrt(u))
 v = (sqrt(3)*((R/2)-y))*rand(1);
 x = sqrt(3)*y + v
% x = 300;
% y = 800;
num = 5;
PPP=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fang算法
%移动台到基站的距离
R1 = distance(X1,Y1,x,y);
R2 = distance(X2,Y2,x,y);
R3 = distance(X3,Y3,x,y);
% R4 = distance(X4,Y4,x,y);
% R5 = distance(X5,Y5,x,y);
% R6 = distance(X6,Y6,x,y);
% R7 = distance(X7,Y7,x,y);

    for i = 1:num
       for k = 1:PPP 
            sgma = i*30;     
            R2_1(k) = Rab(R2,R1) + gngauss(sgma);           
            R3_1(k) = Rab(R3,R1) + gngauss(sgma);           
%             R4_1(k) = Rab(R4,R1) + gngauss(sgma);            
%             R5_1(k) = Rab(R5,R1) + gngauss(sgma);            
%             R6_1(k) = Rab(R6,R1) + gngauss(sgma);            
%             R7_1(k) = Rab(R7,R1) + gngauss(sgma);            
   
  g(k)=(R3_1(k)*X2/R2_1(k)-X3)/Y3;
  h(k)=(K3-R3_1(k)^2+R3_1(k)*R2_1(k)*(1-(X2/R2_1(k))^2))/(2*Y3);
  d(k)=-(1-(X2/R2_1(k))^2+g(k)^2);
  e(k)=X2*(1-(X2/R2_1(k))^2)-2*g(k)*h(k);
  f(k)=(R2_1(k)^2/4)*(1-(X2/R2_1(k))^2)^2-h(k)^2;
%计算得到移动台的位置
  MS_x(k)=(-e(k)-sqrt(e(k)^2-4*d(k)*f(k)))/(2*d(k));
  MS_y(k)=g(k)*MS_x(k)+h(k);
  err(k) = distance(x,y,MS_x(k),MS_y(k));
end
dd(i) = sgma/C;
loc_err3(i)=mean(err);
end

% plot(dd,loc_err7,'bo-'); %7个基站
% hold on;
% plot(dd,loc_err6,'go-'); %6个基站
% hold on;
% plot(dd,loc_err5,'ro-'); %5个基站
% hold on;
% plot(dd,loc_err4,'co-'); %4个基站
% hold on;
plot(dd,loc_err3,'ko-'); %3个基站
xlabel('TDOA误差标准差/s');
ylabel('定位误差均值/m');
axis([0 0.55*10^(-6) 0 500]);