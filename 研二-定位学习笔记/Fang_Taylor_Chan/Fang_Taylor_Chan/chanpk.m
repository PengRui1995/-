%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chan定位算法（也称为两步最小二乘定位算法）
% Chan算法在不同基站参与定位的精度比较
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
C = 300000000; % 光速(米每秒)
R = 1000;      % 小区半径（米）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xi_1 = Xi - X1；Yi_1 = Yi - Y1；Ki = Xi^2 + Yi^2
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
K1 = Ka(X1,Y1);
K2 = Ka(X2,Y2);
K3 = Ka(X3,Y3);
K4 = Ka(X4,Y4);
K5 = Ka(X5,Y5);
K6 = Ka(X6,Y6);
K7 = Ka(X7,Y7);
% 随机产生MS的位置（x，y）
u = rand(1);
x = (R/2)*(1-sqrt(u))
v = (sqrt(3)*((R/2)-x))*rand(1);
y = sqrt(3)*x + v
% x = 300;
% y = 800;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算MS到各基站的距离
R1 = distance(X1,Y1,x,y);
R2 = distance(X2,Y2,x,y);
R3 = distance(X3,Y3,x,y);
R4 = distance(X4,Y4,x,y);
R5 = distance(X5,Y5,x,y);
R6 = distance(X6,Y6,x,y);
R7 = distance(X7,Y7,x,y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPP = 1000;
num1 = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7个基站
    for i = 1:num1
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
            % 假设MS到BSi的距离与MS到BS1的距离差为：Ri_1=Ri-R1,误差的方差为sgma，TDOA（时间差）为：Ti_1＝Ri_1/C
            sgma = i*30;
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
            h = 1/2*[HH(R2_1,X2,Y2,X1,Y1);HH(R3_1,X3,Y3,X1,Y1);HH(R4_1,X4,Y4,X1,Y1);HH(R5_1,X5,Y5,X1,Y1);HH(R6_1,X6,Y6,X1,Y1);HH(R7_1,X7,Y7,X1,Y1)];
            Ga = -[X2_1,Y2_1,Rab(R2,R1);X3_1,Y3_1,Rab(R3,R1);X4_1,Y4_1,Rab(R4,R1);X5_1,Y5_1,Rab(R5,R1);X6_1,Y6_1,Rab(R6,R1);X7_1,Y7_1,Rab(R7,R1)];
            EQ = diag([T2_1_fangcha;T3_1_fangcha;T4_1_fangcha;T5_1_fangcha;T6_1_fangcha;T7_1_fangcha]); 
            Za_1 = inv(Ga'*inv(EQ)*Ga)*Ga'*inv(EQ)*h;                    % Za_1为初始ML估计
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
            Fb = C.^2*B*EQ*B;                                               % Fb是误差矢量，是具有 C.^2*B*EQ*B协方差矩阵的高斯随机矢量
            Za_2 = inv(Ga'*inv(Fb)*Ga)*Ga'*inv(Fb)*h;                       % Za_2为Chan算法的第一次ML估计
            EZa_2 = inv(Ga'*inv(Fb)*Ga);                                    % EZa_2为估计位置Za_2的协方差矩阵
            x00 = Za_2(1);
            y00 = Za_2(2);
            R00 = Za_2(3);
            B_pie = diag([x00 - X1;y00 - Y1;R00]);                           % B'
            Fb_pie = 4*B_pie*EZa_2*B_pie;                                    % Fb'
            h_pie = [(x00 - X1).^2;(y00 - Y1).^2;R00.^2];                    % h' 
            Ga_pie = [1 0;0 1;1 1];                                          % Ga'
            Za_pie = inv(Ga_pie'*inv(Fb_pie)*Ga_pie)*Ga_pie'*inv(Fb_pie)*h_pie;  % Za'
            Zp = sqrt(abs(Za_pie)) + [X1;Y1];                                     % Zp：MS的最终估计位置
            mse(i) = distance(x,y,Zp(1),Zp(2))^2;                            % 求估计位置与真实位置的误差
            tmse(k) = mse(i);
        end
        dd(i) = sgma/C;
        admse_7(i) = sqrt(sum(tmse)/PPP);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6个基站
    for i = 1:num1
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
            % 假设MS到BSi的距离与MS到BS1的距离差为：R2_1=Ri-R1,产生的误差方差为sgma，TDOA（时间差）为：Ti _1＝Ri_1/C
            sgma = i*30;
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
            h = 1/2*[HH(R2_1,X2,Y2,X1,Y1);HH(R3_1,X3,Y3,X1,Y1);HH(R5_1,X5,Y5,X1,Y1);HH(R6_1,X6,Y6,X1,Y1);HH(R7_1,X7,Y7,X1,Y1)];
            Ga = -[X2_1,Y2_1,Rab(R2,R1);X3_1,Y3_1,Rab(R3,R1);X5_1,Y5_1,Rab(R5,R1);X6_1,Y6_1,Rab(R6,R1);X7_1,Y7_1,Rab(R7,R1)];
            EQ = diag([T2_1_fangcha;T3_1_fangcha;T5_1_fangcha;T6_1_fangcha;T7_1_fangcha]);            %EQ为TDOA测量值方差对角阵
%             EQ = cov(Q);                                            % EQ为TDOA协方差矩阵
            Za_1 = inv(Ga'*inv(EQ)*Ga)*Ga'*inv(EQ)*h;                 % Za_1的ML估计
            x0 = Za_1(1);
            y0 = Za_1(2);
            R10 = Za_1(3);
            R20 = distance(X2,Y2,x0,y0);
            R30 = distance(X3,Y3,x0,y0);
%             R40 = distance(X4,Y4,x0,y0);
            R50 = distance(X5,Y5,x0,y0);
            R60 = distance(X6,Y6,x0,y0);
            R70 = distance(X7,Y7,x0,y0);    
            B = diag([R20;R30;R50;R60;R70]);
            Fb = C.^2*B*EQ*B;                                         % Fb是误差矢量，是具有 C.^2*B*EQ*B协方差矩阵的高斯随机矢量
            Za_2 = inv(Ga'*inv(Fb)*Ga)*Ga'*inv(Fb)*h;
            EZa_2 = inv(Ga'*inv(Fb)*Ga);                              % 估计位置Za_2的协方差矩阵
            x00 = Za_2(1);
            y00 = Za_2(2);
            R00 = Za_2(3);
            B_pie = diag([x00 - X1;y00 - Y1;R00]);                     % B'
            Fb_pie = 4*B_pie*EZa_2*B_pie;                              % Fb'
            h_pie = [(x00 - X1).^2;(y00 - Y1).^2;R00.^2];              % h'
            Ga_pie = [1 0;0 1;1 1];                                    % Ga'
            Za_pie = inv(Ga_pie'*inv(Fb_pie)*Ga_pie)*Ga_pie'*inv(Fb_pie)*h_pie;  % Za'
            Zp = sqrt(abs(Za_pie)) + [X1;Y1];                               % Zp：MS的最终估计位置
            mse(i) = distance(x,y,Zp(1),Zp(2))^2;                      % 求估计位置与真实位置的误差
            tmse(k) = mse(i);
        end
        dd(i) = sgma/C;
        admse_6(i) = sqrt(sum(tmse)/PPP);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5个基站
    for i = 1:num1
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
            sgma = i*30;
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
            h = 1/2*[HH(R2_1,X2,Y2,X1,Y1);HH(R3_1,X3,Y3,X1,Y1);HH(R6_1,X6,Y6,X1,Y1);HH(R7_1,X7,Y7,X1,Y1)];
            Ga = -[X2_1,Y2_1,Rab(R2,R1);X3_1,Y3_1,Rab(R3,R1);X6_1,Y6_1,Rab(R6,R1);X7_1,Y7_1,Rab(R7,R1)];
            EQ = diag([T2_1_fangcha;T3_1_fangcha;T6_1_fangcha;T7_1_fangcha]); 
            Za_1 = inv(Ga'*inv(EQ)*Ga)*Ga'*inv(EQ)*h;
            x0 = Za_1(1);
            y0 = Za_1(2);
            R10 = Za_1(3);
            R20 = distance(X2,Y2,x0,y0);
            R30 = distance(X3,Y3,x0,y0);
%             R40 = distance(X4,Y4,x0,y0);
%             R50 = distance(X5,Y5,x0,y0);
            R60 = distance(X6,Y6,x0,y0);
            R70 = distance(X7,Y7,x0,y0);    
            B = diag([R20;R30;R60;R70]);
            Fb = C.^2*B*EQ*B;
            Za_2 = inv(Ga'*inv(Fb)*Ga)*Ga'*inv(Fb)*h;
            EZa_2 = inv(Ga'*inv(Fb)*Ga);
            x00 = Za_2(1);
            y00 = Za_2(2);
            R00 = Za_2(3);
            B_pie = diag([x00 - X1;y00 - Y1;R00]);
            Fb_pie = 4*B_pie*EZa_2*B_pie;
            h_pie = [(x00 - X1).^2;(y00 - Y1).^2;R00.^2];
            Ga_pie = [1 0;0 1;1 1];
            Za_pie = inv(Ga_pie'*inv(Fb_pie)*Ga_pie)*Ga_pie'*inv(Fb_pie)*h_pie;
            Zp = sqrt(abs(Za_pie)) + [X1;Y1];
            mse(i) = distance(x,y,Zp(1),Zp(2))^2;
            tmse(k) = mse(i);
        end
        dd(i) = sgma/C;
        admse_5(i) = sqrt(sum(tmse)/PPP);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4个基站
    for i = 1:num1
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
%             R7_1 = Rab(R7,R1) + gngauss(sgma);
%             T7_1(ss) = R7_1/C;
        end
        T2_1_fangcha = var(T2_1);
        T3_1_fangcha = var(T3_1);
%         T4_1_fangcha = var(T4_1);
%         T5_1_fangcha = var(T5_1);
        T6_1_fangcha = var(T6_1);
%         T7_1_fangcha = var(T7_1);
        for k = 1:PPP
            % 假设MS到BSi的距离与MS到BS1的距离差为：R2_1=Ri-R1,产生的误差方差为sgma，TDOA（时间差）为：Ti_1＝Ri_1/C
            sgma = i*30;
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
%             R7_1 = Rab(R7,R1) + gngauss(sgma);
%             T7_1 = R7_1/C;
            h = 1/2*[HH(R2_1,X2,Y2,X1,Y1);HH(R3_1,X3,Y3,X1,Y1);HH(R6_1,X6,Y6,X1,Y1)];
            Ga = -[X2_1,Y2_1,Rab(R2,R1);X3_1,Y3_1,Rab(R3,R1);X7_1,Y6_1,Rab(R6,R1)];
            EQ = diag([T2_1_fangcha;T3_1_fangcha;T6_1_fangcha]);
            Za_1 = inv(Ga'*inv(EQ)*Ga)*Ga'*inv(EQ)*h;
            x0 = Za_1(1);
            y0 = Za_1(2);
            R10 = Za_1(3);
            R20 = distance(X2,Y2,x0,y0);
            R30 = distance(X3,Y3,x0,y0);
%             R40 = distance(X4,Y4,x0,y0);
%             R50 = distance(X5,Y5,x0,y0);
            R60 = distance(X6,Y6,x0,y0);
%             R70 = distance(X7,Y7,x0,y0);    
            B = diag([R20;R30;R60]);
            Fb = C.^2*B*EQ*B;
            Za_2 = inv(Ga'*inv(Fb)*Ga)*Ga'*inv(Fb)*h;
            EZa_2 = inv(Ga'*inv(Fb)*Ga);
            x00 = Za_2(1);
            y00 = Za_2(2);
            R00 = Za_2(3);
            B_pie = diag([x00 - X1;y00 - Y1;R00]); 
            Fb_pie = 4*B_pie*EZa_2*B_pie; 
            h_pie = [(x00 - X1).^2;(y00 - Y1).^2;R00.^2];
            Ga_pie = [1 0;0 1;1 1];
            Za_pie = inv(Ga_pie'*inv(Fb_pie)*Ga_pie)*Ga_pie'*inv(Fb_pie)*h_pie;
            Zp = sqrt(abs(Za_pie)) + [X1;Y1];
            mse(i) = distance(x,y,Zp(1),Zp(2))^2;
            tmse(k) = mse(i);
        end
        dd(i) = sgma/C;
        admse_4(i) = sqrt(sum(tmse)/PPP);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:num1
        for k = 1:PPP
            % 假设MS到BSi的距离与MS到BS1的距离差为：R2_1=Ri-R1,产生的误差为sgma，TDOA（时间差）为：Ti_1＝Ri_1/C
            sgma = i*30;
            R2_1 = Rab(R2,R1) + gngauss(sgma);
            R7_1 = Rab(R7,R1) + gngauss(sgma);   
            c = [(-R2_1);(-R7_1)];
            d = 1/2*[(K2 - R2_1^2);(K7 - R7_1^2)];
            h = [X2,Y2;X7,Y7];  
            z = R1*inv(h)*c + inv(h)*d  ;
            x_guji(i) = z(1);
            y_guji(i) = z(2);
            mse(i) = distance(x,y,x_guji(i),y_guji(i))^2;
            tmse(k) = mse(i);
        end
        dd(i) = sgma/C;
        admse_3(i) = sqrt(sum(tmse)/PPP);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
plot(dd,admse_7,'bo-'); %7个基站
hold on;
plot(dd,admse_6,'go-'); %6个基站
hold on;
plot(dd,admse_5,'ro-'); %5个基站
hold on;
plot(dd,admse_4,'co-'); %4个基站
hold on;
plot(dd,admse_3,'ko-'); %3个基站
xlabel('TDOA误差标准差/s');
ylabel('定位均方根误差(RMSE)/m');
axis([0 0.55*10^(-6) 0 250]);
legend('7基站','6基站','5基站','4基站','3基站');