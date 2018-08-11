clear
clc
A=[0 1 0 0 0 0 0;
1 0 1 1 0 0 0;
0 1 0 0 0 1 0;
0 1 0 0 1 0 0;
0 0 0 1 0 0 0;
0 0 1 0 0 0 1;
0 0 0 0 0 1 0];
DEF=length(A);
x=100*rand(1,DEF);
y=100*rand(1,DEF);%����ĵ� 
plot(x,y,'ro');      
for i=1:DEF 
    clear a;
    a=find(A(i,:)>0);%��A����ÿ�д���0�������ڸ��еĵ�ַ�ҳ�������a�� 
    if ~isempty(a)
        for j=1:length(a)
            c=num2str(A(i,a(j))); %��A�е�Ȩֵת��Ϊ�ַ���    
            text((x(i)+x(a(j)))/2,(y(i)+y(a(j)))/2,c,'Fontsize',18);%��Ȩֵ��ʾ�����������м�   
            hold on; 
            
            line([x(i) x(a(j))],[y(i) y(a(j))]);%���� 
        end  
    end
end         
title('�������ͼ'); 
e=num2str(DEF); 
legend(e);%���Ͻ���ʾ�ڵ�ĸ��� 
for m=1:DEF 
     A(m,m)=m; 
     f=num2str(A(m,m)); 
     hold on;
     text((x(m)+x(m))/2,(y(m)+y(m))/2,f,'Fontsize',18);%��Ȩֵ��ʾ�����������м�  
end 