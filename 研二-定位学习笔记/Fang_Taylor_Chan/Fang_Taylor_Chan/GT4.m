function y = GT4(R,X,Y,x0,y0);
y = [(X(1) - x0)/R(1) - (X(2) - x0)/R(2),(Y(1) - y0)/R(1) - (Y(2) - y0)/R(2);
     (X(1) - x0)/R(1) - (X(3) - x0)/R(3),(Y(1) - y0)/R(1) - (Y(3) - y0)/R(3);
     (X(1) - x0)/R(1) - (X(4) - x0)/R(4),(Y(1) - y0)/R(1) - (Y(4) - y0)/R(4)];