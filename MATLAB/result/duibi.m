figure(12);
x13=1:1:r1;
y13=1:1:r1;
x23=1:1:r2;
y23=1:1:r2;


for i=1:rmax
    x13(i)=i;
    y13(i) = n - DEAD2(i);
    x23(i)=i;
    y23(i) = n - DEAD22(i);    
end
subplot(2,1,1)
plot(x13,y13,'r--',x23,y23,'b--');
title('A��4������8���������ڵ����Ա�');
xlabel('����������/�֣�');
ylabel('���ڵ�����/����');
legend('8��','4��');


for i=1:rmax
    x13(i)=i;
    y13(i) = STATISTICS.Remain_E2(i);
    x23(i)=i;
    y23(i) = STATISTICS22.Remain_E2(i);;    
end
subplot(2,1,2);
plot(x13,y13,'r--',x23,y23,'b--');
title('B��4������8�����������ĶԱ�');
xlabel('����������/�֣�');
ylabel('ʣ��������/J��');
legend('8��','4��');