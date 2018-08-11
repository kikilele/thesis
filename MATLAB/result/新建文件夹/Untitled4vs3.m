figure(11);
x11=1:1:r1;
y11=1:1:r1;
x21=1:1:r2;
y21=1:1:r2;
x31=1:1:r3;
y31=1:1:r3;


for i=1:rmax
    x11(i)=i;
    y11(i) = n - DEAD2(i);
    x21(i)=i;
    y21(i) = n - DEAD2Copy(i);
    x31(i)=i;
    y31(i) = n - DEAD25(i);
end
subplot(2,1,1);
plot(x11,y11,'r--',x21,y21,'b--',x31,y31,'g--');
title('A£ºSurvival nodes with Closeness Centrality algorithm');
xlabel('Round');
ylabel('Number of survival nodes');
legend('3 sub-regions','4 sub-regions','5 sub-regions');


x31=1:1:r3;
y31=1:1:r3;
x41=1:1:r4;
y41=1:1:r4;
x51=1:1:r5;
y51=1:1:r5;

for i=1:rmax
    
    x31(i)=i;
    y31(i) = n - DEAD4(i); 
    x41(i)=i;
    y41(i) = n - DEAD4Copy(i);
    x51(i)=i;
    y51(i) = n - DEAD45(i);
end
subplot(2,1,2);
plot(x31,y31,'r--',x41,y41,'b--',x51,y51,'g--');
title('A£ºSurvival nodes with PageRank algorithm');
xlabel('Round');
ylabel('Number of survival nodes');
legend('3 sub-regions','4 sub-regions','5 sub-regions');


% figure(11);
% x11=1:1:r1;
% y11=1:1:r1;
% x21=1:1:r2;
% y21=1:1:r2;
% x31=1:1:r3;
% y31=1:1:r3;
% x41=1:1:r4;
% y41=1:1:r4;
% 
% 
% 
% for i=1:rmax
%     x12(i)=i;
%     y12(i) = STATISTICS.Remain_E2(i);
%     x22(i)=i;
%     y22(i) = STATISTICSCopy.Remain_E2(i);
%     x32(i)=i;
%     y32(i) = STATISTICS.Remain_E4(i);
%     x42(i)=i;
%     y42(i) = STATISTICSCopy.Remain_E4(i);
%     
% end
% subplot(2,1,2);
% plot(x12,y12,'r-.',x22,y22,'b-',x32,y32,'g-.',x42,y42,'r-');
% title('B£ºEnergy consumption');
% xlabel('Round');
% ylabel('Residual energy');
% legend('4-cc','3-cc','4-pr','3-pr');