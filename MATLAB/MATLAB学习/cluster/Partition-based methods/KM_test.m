% 调用KM程序实例
clear;
data=rand(100,2);

[best_Label, best_Center, best_ind, label] = KM(data,4,'kmeans');

% [best_Label, best_Center, best_ind, label] = KM(data,4,'kmedoids');

 %画图

for i=1:length(best_Label)
    if (best_Label(i)==1)
        plot(data(i,1),data(i,2),'r.','MarkerSize',8);
        hold on
    elseif (best_Label(i)==2)
        plot(data(i,1),data(i,2),'g.','MarkerSize',8);
        hold on
    elseif (best_Label(i)==3)
        plot(data(i,1),data(i,2),'b.','MarkerSize',8);
        hold on
    else
        plot(data(i,1),data(i,2),'k.','MarkerSize',8);
        hold on
    end
end

for i=1:length(best_Center)
    plot(best_Center(1,i),best_Center(2,i),'kx','MarkerSize',12);
    hold on
end