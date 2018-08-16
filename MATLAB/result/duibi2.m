figure(11);
x11=1:1:r1;
y11=1:1:r1;
x21=1:1:r2;
y21=1:1:r2;
x31=1:1:r3;
y31=1:1:r3;
x51=1:1:r5;
y51=1:1:r5;

for i=1:rmax
    x11(i)=i;
    y11(i) = n - DEAD1(i);
    x21(i)=i;
    y21(i) = n - DEAD2(i);
    x31(i)=i;
    y31(i) = n - DEAD3(i);    
    x51(i)=i;
    y51(i) = n - DEAD5(i);
end
subplot(2,1,1)
plot(x11,y11,'r--',x21,y21,'g--',x31,y31,'b--',x51,y51,'k--');
title('A：存活节点数统计图');
xlabel('工作轮数（/轮）');
ylabel('存活节点数（/个）');
legend('LEACH','NEW','LEACH-C','SEP');

for i=1:rmax
    x12(i)=i;
    y12(i) = STATISTICS.Remain_E1(:,i);
    x22(i)=i;
    y22(i) = STATISTICS.Remain_E2(:,i);
    x32(i)=i;
    y32(i) = STATISTICS.Remain_E3(:,i);
    x52(i)=i;
    y52(i) = STATISTICS.Remain_E5(:,i);
end
subplot(2,1,2)
plot(x12,y12,'r--',x22,y22,'g--',x32,y32,'b--',x52,y52,'k--');
title('B：能力消耗统计图');
xlabel('工作轮数（/轮）');
ylabel('剩余能量（/J）');
legend('LEACH','NEW','LEACH-C','SEP');