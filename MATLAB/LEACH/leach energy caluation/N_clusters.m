el=0.0013*10^-12; es=10*10^-12; Ee=50*10^-9; Ebf=5*10^-9;
n=1000; l=4000; Nf=10000; M=100;
[m]=[1:1:10];
d=150;
f1=el*(d^4)-Ebf.*m-(2*m-1)*Ee;
k=((n/(2*3.14))^(0.5))*((es./f1).^(0.5))*M;
plot(m,k);
xlabel('Head set size', 'FontSize',14)
ylabel('Number Of Clusters', 'FontSize',14)