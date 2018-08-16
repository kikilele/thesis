figure(11);
x11=1:1:r1;
y11=1:1:r1;
x21=1:1:r2;
y21=1:1:r2;


for i=1:rmax
    x11(i)=i;
    y11(i) = n - DEAD24putu(i);
    x21(i)=i;
    y21(i) = n - DEAD24DD(i);
    
end
subplot(2,1,1)
plot(x11,y11,'r--',x21,y21,'b--');
title('A£ºSurvival nodes (4 sub-regions)');
xlabel('Round');
ylabel('Number of survival nodes');
legend('Spectral Graph Theory','Density and Distance');


x12=1:1:r1;
y12=1:1:r1;
x22=1:1:r2;
y22=1:1:r2;


for i=1:rmax
    x12(i)=i;
    y12(i) = n - DEAD25putu(i);
    x22(i)=i;
    y22(i) = n - DEAD25DD(i);
    
end
subplot(2,1,2)
plot(x12,y12,'r--',x22,y22,'b--');
title('A£ºSurvival nodes (5 sub-regions)');
xlabel('Round');
ylabel('Number of survival nodes');
legend('Spectral Graph Theory','Density and Distance');