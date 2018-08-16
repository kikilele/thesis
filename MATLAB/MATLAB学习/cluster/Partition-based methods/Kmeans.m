clc;clear;

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
n=200;


figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    plot(S(i).xd,S(i).yd,'o');
    hold on;
end

X=[XR;YR];
x=X';

% x=[randn(3,2)*.4;randn(4,2)*.5+ones(4,1)*[4 4]];
% x=rand(100,2);
% class1 = kmeans(x, 2);
[class, C, sumd, D] = kmeans(x, 4,'start','cluster');
% class1
class
C

for i=1:length(class)
    if class(i)==1
        plot(x(i,1),x(i,2),'r*');
        plot(C(1,1),C(1,2),'r+');
        hold on;
    elseif class(i)==2
        plot(x(i,1),x(i,2),'b*');
        plot(C(2,1),C(2,2),'b+');
        hold on;
    elseif class(i)==3
        plot(x(i,1),x(i,2),'g*');
        plot(C(3,1),C(3,2),'g+');
        hold on;
    else
        plot(x(i,1),x(i,2),'k*');
        plot(C(4,1),C(4,2),'k+');
        hold on;
    end
end