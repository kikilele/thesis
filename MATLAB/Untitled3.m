clear;

xm=100;
ym=100;
n=10;
figure(1);
for i=1:1:n
    x(i)=rand(1,1)*xm;
    y(i)=rand(1,1)*ym;
    plot(x(i),y(i),'+');
    hold on;      
end

for i=1:n
    for j=1:n
        D(i,j)=sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2 );
        if D(i,j)>50
            D(i,j)=inf;   
        elseif i==j 
            D(i,j)=0;
        else D(i,j)=1;
        end
    end
end

for i=1:n
    clear a;
    a=find(D(i,:)==1);
    if ~isempty(a)
        for j=1:length(a)
            c=num2str(D(i,a(j))); 
            %text((x(i)+x(a(j)))/2,(y(i)+y(a(j)))/2,c,'Fontsize',18);
            hold on;             
            line([x(i) x(a(j))],[y(i) y(a(j))]);
        end  
    end
    e=num2str(i);
    hold on;
    text((x(i)+x(i))/2,(y(i)+y(i))/2,e,'Fontsize',18,'color','red');
end         

M=inf;
p=zeros(length(D));
d=p;
d(:)=M;
for i=1:length(D)    
    pb(1:length(D))=0;
    pb(i)=1;
    index1=i;  
    index2=ones(1,length(D));
    temp=i;    
    d(i,i)=0;
    while sum(pb)<length(D)  %判断是否已经定位完成
    tb=find(pb==0);%查找pb中等于0的元素，即还未定位的那些节点
    d(i,tb)=min(d(i,tb),d(i,temp)+D(temp,tb)); %比较权值大小，将小的存入d
    tmpb=find(d(i,tb)==min(d(i,tb)));%选取最小值的节点进行下一轮的计算
    temp=tb(tmpb(1));%取相同值的第一个值作为下一轮的计算，如果有多个相同值
    pb(temp)=1;
    index1=[index1,temp];%将下一个进行的节点号接入
    index=index1(d(i,index1)==d(i,temp)-D(temp,index1));
    if length(index)>=2
        index=index(1);
    end
    index2(temp)=index;
    end
end
d;



