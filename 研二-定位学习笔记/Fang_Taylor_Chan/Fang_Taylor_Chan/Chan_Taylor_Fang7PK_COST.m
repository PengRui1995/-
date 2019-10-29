%7基站时Chan、Taylor、Fang算法平均定位误差比较
clear;clc;
C = 300000000; % 光速(米每秒)
R = [1000 2000 3000 4000 5000];      % 小区半径（米）
sgma =30;
Rn =1;   %(km)
x_delta = 0;
y_delta = 0;
e = 0;
% 随机产生MS的位置（x，y）
 u = rand(1);
 y = (R/2)*(1-sqrt(u))
 v = (sqrt(3)*((R/2)-y))*rand(1);
 x = sqrt(3)*y + v
num = 5;
PPP=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chan定位算法
    for i = 1:num
% 蜂窝系统各基站坐标
X1 = 0;
Y1 = 0;
X2 = R(i)*(1+sqrt(3)/2);
Y2 = 0;
X3 = R(i)*(sqrt(3)/2);
Y3 = R(i)*(-3/2);
X4 = R(i)*(-sqrt(3)/2);
Y4 = R(i)*(-3/2);
X5 = R(i)*(-1-(sqrt(3)/2));
Y5 = 0;
X6 = R(i)*(-sqrt(3)/2);
Y6 = R(i)*(3/2);
X7 = R(i)*(sqrt(3)/2);
Y7 = R(i)*(3/2);
% 计算MS到各基站的距离
R1 = distance(X1,Y1,x(i),y(i));
R2 = distance(X2,Y2,x(i),y(i));
R3 = distance(X3,Y3,x(i),y(i));
R4 = distance(X4,Y4,x(i),y(i));
R5 = distance(X5,Y5,x(i),y(i));
R6 = distance(X6,Y6,x(i),y(i));
R7 = distance(X7,Y7,x(i),y(i));
 

% KM = XM^2 + YM^2
K1 = Ka(X1,Y1);
K2 = Ka(X2,Y2);
K3 = Ka(X3,Y3);
K4 = Ka(X4,Y4);
K5 = Ka(X5,Y5);
K6 = Ka(X6,Y6);
K7 = Ka(X7,Y7);
% XM_1 = XM - X1；YM_1 = YM - X1
X2_1 = Xab(X2,X1);
Y2_1 = Xab(Y2,Y1);
X3_1 = Xab(X3,X1);
Y3_1 = Xab(Y3,Y1);
X4_1 = Xab(X4,X1);
Y4_1 = Xab(Y4,Y1);
X5_1 = Xab(X5,X1);
Y5_1 = Xab(Y5,Y1);
X6_1 = Xab(X6,X1);
Y6_1 = Xab(Y6,Y1);
X7_1 = Xab(X7,X1);
Y7_1 = Xab(Y7,Y1);
 t_sgma=10^(-7);   
 % 各TDOA测量值的方差
            D_med1= 1.82*(0.4*10^(-6)*(R1/1000)^0.5)^2;
            D_med2= 1.82*(0.4*10^(-6)*(R2/1000)^0.5)^2;
            D_med3= 1.82*(0.4*10^(-6)*(R3/1000)^0.5)^2;
            D_med4= 1.82*(0.4*10^(-6)*(R4/1000)^0.5)^2;
            D_med5= 1.82*(0.4*10^(-6)*(R5/1000)^0.5)^2;
            D_med6= 1.82*(0.4*10^(-6)*(R6/1000)^0.5)^2;
            D_med7= 1.82*(0.4*10^(-6)*(R7/1000)^0.5)^2;
        T2_1_fangcha = D_med2+D_med1+2*t_sgma^2;
        T3_1_fangcha = D_med3+D_med1+2*t_sgma^2;
        T4_1_fangcha = D_med4+D_med1+2*t_sgma^2;
        T5_1_fangcha = D_med5+D_med1+2*t_sgma^2;
        T6_1_fangcha = D_med6+D_med1+2*t_sgma^2;
        T7_1_fangcha = D_med7+D_med1+2*t_sgma^2;
        
        for k = 1:PPP
            t_rms1=0.4*10^(-6)*(R1/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms2=0.4*10^(-6)*(R2/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms3=0.4*10^(-6)*(R3/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms4=0.4*10^(-6)*(R4/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms5=0.4*10^(-6)*(R5/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms6=0.4*10^(-6)*(R6/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            t_rms7=0.4*10^(-6)*(R7/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            % 假设MS到BSi的距离与MS到BS1的距离差为：R2_1=Ri-R1,产生的误差方差为sgma，TDOA（时间差）为：Ti_1＝Ri_1/C
            R2_1 = Rab(R2,R1) + gngauss(sgma)+C*(t_rms2-t_rms1);
            R3_1 = Rab(R3,R1) + gngauss(sgma)+C*(t_rms3-t_rms1);
            R4_1 = Rab(R4,R1) + gngauss(sgma)+C*(t_rms4-t_rms1);
            R5_1 = Rab(R5,R1) + gngauss(sgma)+C*(t_rms5-t_rms1);
            R6_1 = Rab(R6,R1) + gngauss(sgma)+C*(t_rms6-t_rms1);
            R7_1 = Rab(R7,R1) + gngauss(sgma)+C*(t_rms7-t_rms1);

            h = 1/2*[HH(R2_1,X2,Y2,X1,Y1);HH(R3_1,X3,Y3,X1,Y1);HH(R4_1,X4,Y4,X1,Y1);HH(R5_1,X5,Y5,X1,Y1);HH(R6_1,X6,Y6,X1,Y1);HH(R7_1,X7,Y7,X1,Y1)];
            Ga = -[X2_1,Y2_1,R2_1;X3_1,Y3_1,R3_1;X4_1,Y4_1,R4_1;X5_1,Y5_1,R5_1;X6_1,Y6_1,R6_1;X7_1,Y7_1,R7_1];
            EQ = diag([T2_1_fangcha;T3_1_fangcha;T4_1_fangcha;T5_1_fangcha;T6_1_fangcha;T7_1_fangcha]);  % Q矩阵为一对角阵
            % Za的第一次WLS估计
            Za_1 = inv(Ga'*inv(EQ)*Ga)*Ga'*inv(EQ)*h; 
            x0 = Za_1(1);
            y0 = Za_1(2);
            R10 = Za_1(3);
            R20 = distance(X2,Y2,x0,y0);
            R30 = distance(X3,Y3,x0,y0);
            R40 = distance(X4,Y4,x0,y0);
            R50 = distance(X5,Y5,x0,y0);
            R60 = distance(X6,Y6,x0,y0);
            R70 = distance(X7,Y7,x0,y0);    
            B = diag([R20;R30;R40;R50;R60;R70]); 
            Fb = C.^2*B*EQ*B;
            % Za的第二次WLS估计
            Za_2 = inv(Ga'*inv(Fb)*Ga)*Ga'*inv(Fb)*h;
            Ga0 = -[X2_1,Y2_1,Rab(R20,R10);X3_1,Y3_1,Rab(R30,R10);X4_1,Y4_1,Rab(R40,R10);X5_1,Y5_1,Rab(R50,R10);X6_1,Y6_1,Rab(R60,R10);X7_1,Y7_1,Rab(R70,R10)];
            EZa_2 = inv(Ga0'*inv(Fb)*Ga0); 
            x00 = Za_2(1);
            y00 = Za_2(2);
            R00 = Za_2(3);
            B_pie = diag([x00 - X1;y00 - Y1;R00]);
            Fb_pie = 4*B_pie*EZa_2*B_pie;
            h_pie = [(x00 - X1).^2;(y00 - Y1).^2;R00.^2]; 
            Ga_pie = [1 0;0 1;1 1];
            % Za的第三次WLS估计
            Za_pie = inv(Ga_pie'*inv(Fb_pie)*Ga_pie)*Ga_pie'*inv(Fb_pie)*h_pie;
            Zp = sqrt(abs(Za_pie)) + [X1;Y1];
        err_chan(k) = distance(x(i),y(i),Zp(1),Zp(2));
        end
        loc_err_Chan7(i) = mean(err_chan);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taylor定位算法
for i = 1:num
% 蜂窝系统各基站坐标
X1 = 0;
Y1 = 0;
X2 = R(i)*(1+sqrt(3)/2);
Y2 = 0;
X3 = R(i)*(sqrt(3)/2);
Y3 = R(i)*(-3/2);
X4 = R(i)*(-sqrt(3)/2);
Y4 = R(i)*(-3/2);
X5 = R(i)*(-1-(sqrt(3)/2));
Y5 = 0;
X6 = R(i)*(-sqrt(3)/2);
Y6 = R(i)*(3/2);
X7 = R(i)*(sqrt(3)/2);
Y7 = R(i)*(3/2);
% 计算MS到各基站的距离
R1 = distance(X1,Y1,x(i),y(i));
R2 = distance(X2,Y2,x(i),y(i));
R3 = distance(X3,Y3,x(i),y(i));
R4 = distance(X4,Y4,x(i),y(i));
R5 = distance(X5,Y5,x(i),y(i));
R6 = distance(X6,Y6,x(i),y(i));
R7 = distance(X7,Y7,x(i),y(i));
 
% KM = XM^2 + YM^2
K1 = Ka(X1,Y1);
K2 = Ka(X2,Y2);
K3 = Ka(X3,Y3);
K4 = Ka(X4,Y4);
K5 = Ka(X5,Y5);
K6 = Ka(X6,Y6);
K7 = Ka(X7,Y7);
% XM_1 = XM - X1；YM_1 = YM - X1
X2_1 = Xab(X2,X1);
Y2_1 = Xab(Y2,Y1);
X3_1 = Xab(X3,X1);
Y3_1 = Xab(Y3,Y1);
X4_1 = Xab(X4,X1);
Y4_1 = Xab(Y4,Y1);
X5_1 = Xab(X5,X1);
Y5_1 = Xab(Y5,Y1);
X6_1 = Xab(X6,X1);
Y6_1 = Xab(Y6,Y1);
X7_1 = Xab(X7,X1);
Y7_1 = Xab(Y7,Y1);
 t_sgma=10^(-7);   
 % 各TDOA测量值的方差
            D_med1= 1.82*(0.4*10^(-6)*(R1/1000)^0.5)^2;
            D_med2= 1.82*(0.4*10^(-6)*(R2/1000)^0.5)^2;
            D_med3= 1.82*(0.4*10^(-6)*(R3/1000)^0.5)^2;
            D_med4= 1.82*(0.4*10^(-6)*(R4/1000)^0.5)^2;
            D_med5= 1.82*(0.4*10^(-6)*(R5/1000)^0.5)^2;
            D_med6= 1.82*(0.4*10^(-6)*(R6/1000)^0.5)^2;
            D_med7= 1.82*(0.4*10^(-6)*(R7/1000)^0.5)^2;
        T2_1_fangcha = D_med2+D_med1+2*t_sgma^2;
        T3_1_fangcha = D_med3+D_med1+2*t_sgma^2;
        T4_1_fangcha = D_med4+D_med1+2*t_sgma^2;
        T5_1_fangcha = D_med5+D_med1+2*t_sgma^2;
        T6_1_fangcha = D_med6+D_med1+2*t_sgma^2;
        T7_1_fangcha = D_med7+D_med1+2*t_sgma^2;
        for k = 1:PPP
            t_rms1=0.4*10^(-6)*(R1/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms2=0.4*10^(-6)*(R2/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms3=0.4*10^(-6)*(R3/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms4=0.4*10^(-6)*(R4/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms5=0.4*10^(-6)*(R5/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms6=0.4*10^(-6)*(R6/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms7=0.4*10^(-6)*(R7/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 
            % 假设初始位置为x0，y0
            x0 = x(i) + 10;%gngauss(sgma);
            y0 = y(i) + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma)+C*(t_rms2-t_rms1);
            T2_1 = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma)+C*(t_rms3-t_rms1);
            T3_1 = R3_1/C;
            R4_1 = Rab(R4,R1) + gngauss(sgma)+C*(t_rms4-t_rms1);
            T4_1 = R4_1/C;
            R5_1 = Rab(R5,R1) + gngauss(sgma)+C*(t_rms5-t_rms1);
            T5_1 = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma)+C*(t_rms6-t_rms1);
            T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma)+C*(t_rms7-t_rms1);
            T7_1 = R7_1/C;
            if (-Rn*100 < x0 < R(i)/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R(i)/2 + Rn*100) % 设置一个收敛范围
                while abs(e) < R(i)^2
                % MS的估计位置到BSi之间的距离
                R10 = distance(X1,Y1,x0,y0);
                R20 = distance(X2,Y2,x0,y0);
                R30 = distance(X3,Y3,x0,y0);
                R40 = distance(X4,Y4,x0,y0);
                R50 = distance(X5,Y5,x0,y0);
                R60 = distance(X6,Y6,x0,y0);
                R70 = distance(X7,Y7,x0,y0);

                ht =  [R2_1 - (R20 - R10);
                       R3_1 - (R30 - R10);
                       R4_1 - (R40 - R10);
                       R5_1 - (R50 - R10);
                       R6_1 - (R60 - R10);
                       R7_1 - (R70 - R10)];
                EQ = diag([T2_1_fangcha;T3_1_fangcha;T4_1_fangcha;T5_1_fangcha;T6_1_fangcha;T7_1_fangcha]);
                RR = [R10,R20,R30,R40,R50,R60,R70];
                XX = [X1,X2,X3,X4,X5,X6,X7];
                YY = [Y1,Y2,Y3,Y4,Y5,Y6,Y7];
                Gt = GT7(RR,XX,YY,x0,y0);
                delta = inv(Gt'*inv(EQ)*Gt)*Gt'*inv(EQ)*ht;
                x_delta = delta(1);
                y_delta = delta(2);
                x0 = x0 + x_delta;
                y0 = y0 + y_delta;
                e = abs(x_delta) + abs(y_delta);

                    if e < 1;
                        break;
                   end
               end
           end
            err_taylor(k) = distance(x(i),y(i),x0,y0);
        end
        loc_err_Taylor7(i) = mean(err_taylor);   
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Fang定位算法 
for i = 1:num
% 蜂窝系统各基站坐标
X1 = 0;
Y1 = 0;
X2 = R(i)*(1+sqrt(3)/2);
Y2 = 0;
X3 = R(i)*(sqrt(3)/2);
Y3 = R(i)*(-3/2);
X4 = R(i)*(-sqrt(3)/2);
Y4 = R(i)*(-3/2);
X5 = R(i)*(-1-(sqrt(3)/2));
Y5 = 0;
X6 = R(i)*(-sqrt(3)/2);
Y6 = R(i)*(3/2);
X7 = R(i)*(sqrt(3)/2);
Y7 = R(i)*(3/2);
% 计算MS到各基站的距离
R1 = distance(X1,Y1,x(i),y(i));
R2 = distance(X2,Y2,x(i),y(i));
R3 = distance(X3,Y3,x(i),y(i));
R4 = distance(X4,Y4,x(i),y(i));
R5 = distance(X5,Y5,x(i),y(i));
R6 = distance(X6,Y6,x(i),y(i));
R7 = distance(X7,Y7,x(i),y(i));


% KM = XM^2 + YM^2
K1 = Ka(X1,Y1);
K2 = Ka(X2,Y2);
K3 = Ka(X3,Y3);
K4 = Ka(X4,Y4);
K5 = Ka(X5,Y5);
K6 = Ka(X6,Y6);
K7 = Ka(X7,Y7);
% XM_1 = XM - X1；YM_1 = YM - X1
X2_1 = Xab(X2,X1);
Y2_1 = Xab(Y2,Y1);
X3_1 = Xab(X3,X1);
Y3_1 = Xab(Y3,Y1);
X4_1 = Xab(X4,X1);
Y4_1 = Xab(Y4,Y1);
X5_1 = Xab(X5,X1);
Y5_1 = Xab(Y5,Y1);
X6_1 = Xab(X6,X1);
Y6_1 = Xab(Y6,Y1);
X7_1 = Xab(X7,X1);
Y7_1 = Xab(Y7,Y1); 
       for k = 1:PPP  
           t_rms1=0.4*10^(-6)*(R1/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms2=0.4*10^(-6)*(R2/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
 t_rms3=0.4*10^(-6)*(R3/1000)^0.5*exp(gngauss(0,4)/(10*log10(10)));
            R2_1(k) = Rab(R2,R1) + gngauss(sgma)+C*(t_rms2-t_rms1);           
            R3_1(k) = Rab(R3,R1) + gngauss(sgma)+C*(t_rms3-t_rms1);           
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
  err_fang(k) = distance(x(i),y(i),MS_x(k),MS_y(k));
end
loc_err_Fang7(i)=mean(err_fang);
end
%7个基站
plot(R,loc_err_Chan7,'k+-'); hold on; 
plot(R,loc_err_Taylor7,'kd-','MarkerFaceColor','k'); hold on; 
plot(R,loc_err_Fang7,'k^-','MarkerFaceColor','k');
xlabel('基站小区半径');
ylabel('定位误差均值/m');
axis([0 5500 0 500]);
legend('Chan','Taylor','Fang');