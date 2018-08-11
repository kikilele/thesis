el=0.0013*10^-12; es=10*10^-12; Ee=50*10^-9; Ebf=5*10^-9;
n=1000; k=50; d=100; l=4000; r=250*10^3; Estart=0.1;
[M,m]=meshgrid([1:10:200],[1:1:20]);
Ech_elect=l.*(Ee*(n/k)+Ebf*((n/k)-1)+es*d^2);
Enon_ch_elect=l.*(Ee*(1+k)+(k*Ebf)+es.*((M.*M)/(2*3.14159*k)));
Ech_frame=l.*(((n/k)-m)*Ebf+(el*d^4)+((n/k)-m+1)*Ee);
Enon_ch_frame=l*(Ee+es.*((M.^2)/(2*3.14159*k)));
f1=1./(((n/k)-m+1)*k); f2=((n/k)-m)./(((n/k)-m+1)*k);
p1=((m.*Estart)-Ech_elect-Enon_ch_elect);
p2=(f1.*Ech_frame)+(f2.*Enon_ch_frame);
Nf=p1./p2;
surf(M,m,Nf);
xlabel('Network Diameter (m)', 'FontSize',14)
ylabel('Head-Set Size', 'FontSize',14)
zlabel('# of Frames','FontSize',14)