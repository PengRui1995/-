%不同基站时Taylor算法定位性能比较
clear;clc;
C = 300000000; % 光速(米每秒)
R = 1000;      % 小区半径（米）
Rn =1;   %(km)
x_delta = 0;
y_delta = 0;
e = 0;
% stda = 0.01;
% 蜂窝系统各基站坐标
X1 = 0;
Y1 = 0;
X2 = 0;
Y2 = R*(1+sqrt(3)/2);
X3 = R*(-3/2);
Y3 = R*(sqrt(3)/2);
X4 = R*(-3/2);
Y4 = R*(-sqrt(3)/2);
X5 = 0;
Y5 = R*(-1+(sqrt(3)/2));
X6 = R*(3/2);
Y6 = R*(-sqrt(3)/2);
X7 = R*(3/2);
Y7 = R*(sqrt(3)/2);

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
% 随机产生MS的位置（x，y）
 u = rand(1);
 x = (R/2)*(1-sqrt(u))
 v = (sqrt(3)*((R/2)-x))*rand(1);
 y = sqrt(3)*x + v
% x = 300
% y = 800
% 计算MS到各基站的距离
R1 = distance(X1,Y1,x,y);
R2 = distance(X2,Y2,x,y);
R3 = distance(X3,Y3,x,y);
R4 = distance(X4,Y4,x,y);
R5 = distance(X5,Y5,x,y);
R6 = distance(X6,Y6,x,y);
R7 = distance(X7,Y7,x,y);

num = 5;
PPP=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7个基站
for i = 1:num
        sgma = i*30;
        for ss = 1:PPP
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1(ss) = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1(ss) = R3_1/C;
            R4_1 = Rab(R4,R1) + gngauss(sgma);
            T4_1(ss) = R4_1/C;
            R5_1 = Rab(R5,R1) + gngauss(sgma);
            T5_1(ss) = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1(ss) = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1(ss) = R7_1/C;
        end
        
        T2_1_fangcha = var(T2_1);
        T3_1_fangcha = var(T3_1);
        T4_1_fangcha = var(T4_1);
        T5_1_fangcha = var(T5_1);
        T6_1_fangcha = var(T6_1);
        T7_1_fangcha = var(T7_1);
        
        for k = 1:PPP
            sgma = 30*i;
            % 假设初始位置为x0，y0
            x0 = x + 10;%gngauss(sgma);
            y0 = y + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1 = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1 = R3_1/C;
            R4_1 = Rab(R4,R1) + gngauss(sgma);
            T4_1 = R4_1/C;
            R5_1 = Rab(R5,R1) + gngauss(sgma);
            T5_1 = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1 = R7_1/C;

            if (-Rn*100 < x0 < R/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R/2 + Rn*100) % 设置一个收敛范围

                while abs(e) < R.^2
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
            err_taylor(k) = distance(x,y,x0,y0)^2;
        end
        dd(i) = sgma/C;
        admse_taylor7(i) = sqrt(sum(err_taylor)/PPP);   
   end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %6个基站
 for i = 1:num
        sgma = i*30;
        for ss = 1:PPP
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1(ss) = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1(ss) = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1(ss) = R4_1/C;
            R5_1 = Rab(R5,R1) + gngauss(sgma);
            T5_1(ss) = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1(ss) = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1(ss) = R7_1/C;
        end
        
        T2_1_fangcha = var(T2_1);
        T3_1_fangcha = var(T3_1);
%         T4_1_fangcha = var(T4_1);
        T5_1_fangcha = var(T5_1);
        T6_1_fangcha = var(T6_1);
        T7_1_fangcha = var(T7_1);
        
        for k = 1:PPP
            sgma = 30*i;
            % 假设初始位置为x0，y0
            x0 = x + 10;%gngauss(sgma);
            y0 = y + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1 = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1 = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1 = R4_1/C;
            R5_1 = Rab(R5,R1) + gngauss(sgma);
            T5_1 = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1 = R7_1/C;

            if (-Rn*100 < x0 < R/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R/2 + Rn*100) % 设置一个收敛范围

                while abs(e) < R.^2
                % MS的估计位置到BSi之间的距离
                R10 = distance(X1,Y1,x0,y0);
                R20 = distance(X2,Y2,x0,y0);
                R30 = distance(X3,Y3,x0,y0);
%                 R40 = distance(X4,Y4,x0,y0);
                R50 = distance(X5,Y5,x0,y0);
                R60 = distance(X6,Y6,x0,y0);
                R70 = distance(X7,Y7,x0,y0);

                ht =  [R2_1 - (R20 - R10);
                       R3_1 - (R30 - R10);
                       R5_1 - (R50 - R10);
                       R6_1 - (R60 - R10);
                       R7_1 - (R70 - R10)];
                EQ = diag([T2_1_fangcha;T3_1_fangcha;T5_1_fangcha;T6_1_fangcha;T7_1_fangcha]);
                RR = [R10,R20,R30,R50,R60,R70];
                XX = [X1,X2,X3,X5,X6,X7];
                YY = [Y1,Y2,Y3,Y5,Y6,Y7];
                Gt = GT6(RR,XX,YY,x0,y0);
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
           err_taylor(k) = distance(x,y,x0,y0)^2;
        end
        dd(i) = sgma/C;
        admse_taylor6(i) = sqrt(sum(err_taylor)/PPP);   
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %5个基站
 for i = 1:num
        sgma = i*30;
        for ss = 1:PPP
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1(ss) = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1(ss) = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1(ss) = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1(ss) = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1(ss) = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1(ss) = R7_1/C;
        end
        
        T2_1_fangcha = var(T2_1);
        T3_1_fangcha = var(T3_1);
%         T4_1_fangcha = var(T4_1);
%         T5_1_fangcha = var(T5_1);
        T6_1_fangcha = var(T6_1);
        T7_1_fangcha = var(T7_1);
        
        for k = 1:PPP
            sgma = 30*i;
            % 假设初始位置为x0，y0
            x0 = x + 10;%gngauss(sgma);
            y0 = y + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1 = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1 = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1 = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1 = R5_1/C;
            R6_1 = Rab(R6,R1) + gngauss(sgma);
            T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1 = R7_1/C;

            if (-Rn*100 < x0 < R/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R/2 + Rn*100) % 设置一个收敛范围

                while abs(e) < R.^2
                % MS的估计位置到BSi之间的距离
                R10 = distance(X1,Y1,x0,y0);
                R20 = distance(X2,Y2,x0,y0);
                R30 = distance(X3,Y3,x0,y0);
%                 R40 = distance(X4,Y4,x0,y0);
%                 R50 = distance(X5,Y5,x0,y0);
                R60 = distance(X6,Y6,x0,y0);
                R70 = distance(X7,Y7,x0,y0);

                ht =  [R2_1 - (R20 - R10);
                       R3_1 - (R30 - R10);
                       R6_1 - (R60 - R10);
                       R7_1 - (R70 - R10)];
                EQ = diag([T2_1_fangcha;T3_1_fangcha;T6_1_fangcha;T7_1_fangcha]);
                RR = [R10,R20,R30,R60,R70];
                XX = [X1,X2,X3,X6,X7];
                YY = [Y1,Y2,Y3,Y6,Y7];
                Gt = GT5(RR,XX,YY,x0,y0);
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
           err_taylor(k) = distance(x,y,x0,y0)^2;
        end
        dd(i) = sgma/C;
        admse_taylor5(i) = sqrt(sum(err_taylor)/PPP);   
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %4个基站
 for i = 1:num
        sgma = i*30;
        for ss = 1:PPP
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1(ss) = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1(ss) = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1(ss) = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1(ss) = R5_1/C;
%             R6_1 = Rab(R6,R1) + gngauss(sgma);
%             T6_1(ss) = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1(ss) = R7_1/C;
        end
        
        T2_1_fangcha = var(T2_1);
        T3_1_fangcha = var(T3_1);
%         T4_1_fangcha = var(T4_1);
%         T5_1_fangcha = var(T5_1);
%         T6_1_fangcha = var(T6_1);
        T7_1_fangcha = var(T7_1);
        
        for k = 1:PPP
            sgma = 30*i;
            % 假设初始位置为x0，y0
            x0 = x + 10;%gngauss(sgma);
            y0 = y + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1 = R2_1/C;
            R3_1 = Rab(R3,R1) + gngauss(sgma);
            T3_1 = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1 = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1 = R5_1/C;
%             R6_1 = Rab(R6,R1) + gngauss(sgma);
%             T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1 = R7_1/C;

            if (-Rn*100 < x0 < R/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R/2 + Rn*100) % 设置一个收敛范围

                while abs(e) < R.^2
                % MS的估计位置到BSi之间的距离
                R10 = distance(X1,Y1,x0,y0);
                R20 = distance(X2,Y2,x0,y0);
                R30 = distance(X3,Y3,x0,y0);
%                 R40 = distance(X4,Y4,x0,y0);
%                 R50 = distance(X5,Y5,x0,y0);
%                 R60 = distance(X6,Y6,x0,y0);
                R70 = distance(X7,Y7,x0,y0);

                ht =  [R2_1 - (R20 - R10);
                       R3_1 - (R30 - R10);
                       R7_1 - (R70 - R10)];
                EQ = diag([T2_1_fangcha;T3_1_fangcha;T7_1_fangcha]);
                RR = [R10,R20,R30,R70];
                XX = [X1,X2,X3,X7];
                YY = [Y1,Y2,Y3,Y7];
                Gt = GT4(RR,XX,YY,x0,y0);
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
           err_taylor(k) = distance(x,y,x0,y0)^2;
        end
        dd(i) = sgma/C;
        admse_taylor4(i) = sqrt(sum(err_taylor)/PPP);   
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %3个基站
 for i = 1:num
        sgma = i*30;
        for ss = 1:PPP
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1(ss) = R2_1/C;
%             R3_1 = Rab(R3,R1) + gngauss(sgma);
%             T3_1(ss) = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1(ss) = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1(ss) = R5_1/C;
%             R6_1 = Rab(R6,R1) + gngauss(sgma);
%             T6_1(ss) = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1(ss) = R7_1/C;
        end
        
        T2_1_fangcha = var(T2_1);
%         T3_1_fangcha = var(T3_1);
%         T4_1_fangcha = var(T4_1);
%         T5_1_fangcha = var(T5_1);
%         T6_1_fangcha = var(T6_1);
        T7_1_fangcha = var(T7_1);
        
        for k = 1:PPP
            sgma = 30*i;
            % 假设初始位置为x0，y0
            x0 = x + 10;%gngauss(sgma);
            y0 = y + 10;%gngauss(sgma);

            R2_1 = Rab(R2,R1) + gngauss(sgma);
            T2_1 = R2_1/C;
%             R3_1 = Rab(R3,R1) + gngauss(sgma);
%             T3_1 = R3_1/C;
%             R4_1 = Rab(R4,R1) + gngauss(sgma);
%             T4_1 = R4_1/C;
%             R5_1 = Rab(R5,R1) + gngauss(sgma);
%             T5_1 = R5_1/C;
%             R6_1 = Rab(R6,R1) + gngauss(sgma);
%             T6_1 = R6_1/C;
            R7_1 = Rab(R7,R1) + gngauss(sgma);
            T7_1 = R7_1/C;

            if (-Rn*100 < x0 < R/2 + Rn*100) && (-Rn*100 < y0 < sqrt(3)*R/2 + Rn*100) % 设置一个收敛范围

                while abs(e) < R.^2
                % MS的估计位置到BSi之间的距离
                R10 = distance(X1,Y1,x0,y0);
                R20 = distance(X2,Y2,x0,y0);
%                 R30 = distance(X3,Y3,x0,y0);
%                 R40 = distance(X4,Y4,x0,y0);
%                 R50 = distance(X5,Y5,x0,y0);
%                 R60 = distance(X6,Y6,x0,y0);
                R70 = distance(X7,Y7,x0,y0);

                ht =  [R2_1 - (R20 - R10);
                       R7_1 - (R70 - R10)];
                EQ = diag([T2_1_fangcha;T7_1_fangcha]);
                RR = [R10,R20,R70];
                XX = [X1,X2,X7];
                YY = [Y1,Y2,Y7];
                Gt = GT3(RR,XX,YY,x0,y0);
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
           err_taylor(k) = distance(x,y,x0,y0)^2;
        end
        dd(i) = sgma/C;
        admse_taylor3(i) = sqrt(sum(err_taylor)/PPP);   
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(dd,admse_taylor7,'bo-'); %7个基站
hold on;
plot(dd,admse_taylor6,'go-'); %6个基站
hold on;
plot(dd,admse_taylor5,'ro-'); %5个基站
hold on;
plot(dd,admse_taylor4,'co-'); %4个基站
hold on;
plot(dd,admse_taylor3,'ko-'); %3个基站
xlabel('TDOA误差标准差/s');
ylabel('定位均方根误差(RMSE)/m');
axis([0 0.55*10^(-6) 0 250]);
legend('7基站','6基站','5基站','4基站','3基站');