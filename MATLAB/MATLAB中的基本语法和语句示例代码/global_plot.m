function global_plot( )
global X         %����XΪȫ�ֱ���
X=0:0.1:2*pi;
plot_sin(2)
plot_cos(2)


function plot_sin(a)
global X         %ʹ��ȫ�ֱ���XʱҲҪ��global����
y=a*sin(X);
figure
plot(X,y)

function plot_cos(a)
global X         %ʹ��ȫ�ֱ���XʱҲҪ��global����
X=-pi:0.1:pi;    %ȫ�ֱ������޸�
y=a*cos(X);
figure
plot(X,y)
