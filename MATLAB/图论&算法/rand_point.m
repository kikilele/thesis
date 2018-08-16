
clear;
n=10;
xm=30;

figure(1);
for i=1:1:n
    XR(i,1)=rand(1,1)*xm;    
    XR(i,2)=rand(1,1)*xm;
    plot(XR(i,1),XR(i,2),'o');
    hold on;
end

linjiejuzheng=pdist2(XR,XR);

for i=1:length(linjiejuzheng)
    for j=1:length(linjiejuzheng)
        if (linjiejuzheng(i,j) > 10)
            linjiejuzheng(i,j)=inf;
        else
            linjiejuzheng(i,j)=1;
        end
    end
end

D = shortdf(linjiejuzheng);

for i=1:n
    for j=1:n
        if linjiejuzheng(i,j) ~= inf 
            c=num2str(linjiejuzheng(i,j));
            line([XR(i,1) XR(j,1)],[XR(i,2) XR(j,2)]);
        end
        text(XR(i,1),XR(i,2),num2str(i),'Fontsize',14,'color','r');
        hold on
    end
end
