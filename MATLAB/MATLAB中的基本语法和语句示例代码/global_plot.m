function global_plot( )
global X         %定义X为全局变量
X=0:0.1:2*pi;
plot_sin(2)
plot_cos(2)


function plot_sin(a)
global X         %使用全局变量X时也要用global定义
y=a*sin(X);
figure
plot(X,y)

function plot_cos(a)
global X         %使用全局变量X时也要用global定义
X=-pi:0.1:pi;    %全局变量被修改
y=a*cos(X);
figure
plot(X,y)
