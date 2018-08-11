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
y=100*rand(1,DEF);%随机的点 
plot(x,y,'ro');      
for i=1:DEF 
    clear a;
    a=find(A(i,:)>0);%将A矩阵每行大于0的数的在该行的地址找出来放在a中 
    if ~isempty(a)
        for j=1:length(a)
            c=num2str(A(i,a(j))); %将A中的权值转化为字符型    
            text((x(i)+x(a(j)))/2,(y(i)+y(a(j)))/2,c,'Fontsize',18);%将权值显示在两点连线中间   
            hold on; 
            
            line([x(i) x(a(j))],[y(i) y(a(j))]);%连线 
        end  
    end
end         
title('随机拓扑图'); 
e=num2str(DEF); 
legend(e);%左上角显示节点的个数 
for m=1:DEF 
     A(m,m)=m; 
     f=num2str(A(m,m)); 
     hold on;
     text((x(m)+x(m))/2,(y(m)+y(m))/2,f,'Fontsize',18);%将权值显示在两点连线中间  
end 