%����������Բ�������������������ͼ�Լ���Ӧ���ڽӾ���
%����Ϊnum������

function [L]=connectionM(num)
%��������A������ű�ʾ����֮��Ȩֵ�ľ���A��Ҳ�����ڽӾ�����ô����֮��Ȩֵ��Ϊ��Ԫ�صĸ�����Ϊ�õ�Ķ���
DEF=num;
L=rand(DEF);%�������һ��DEF�׾���Ȼ���Դ�Ϊ�����޸�Ϊ�ڽӾ���
for i=1:DEF
    L(i,i)=0;%�Խ�����Ԫ�ظ�Ϊ0������ͼ���ڽӾ��������
end
L=10*L;
L=floor(L);%��ȡ��
L=mod(L,2);%�ڽӾ����е�Ԫ��ȡֵΪ0��1��2
for i=1:DEF
    for j=1:i
        L(j,i)=L(i,j);
    end
end
x=100*rand(1,DEF);
y=100*rand(1,DEF);
plot(x,y,'r+');
for i=1:DEF
    a=find(L(i,:)>0);%��L����ÿ�д���0�������ڸ��еĵ�ַ�ҳ�������a��
    for j=1:length(a)        
        hold on;
        line([x(i),x(a(j))],[y(i),y(a(j))]);
    end
end
title('�������ͼ');
e=num2str(DEF); 
legend(e);
for m=1:DEF
    f=num2str(m); 
    hold on;
    text((x(m)+x(m))/2,(y(m)+y(m))/2,f,'Fontsize',18); 
end
hold off