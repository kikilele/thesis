el=0.0013*10^-12; es=10*10^-12; Ee=50*10^-9; Ebf=5*10^-9;
n=1000; k=50; d=100; l=4000; Nf=10000;
[M,m]=meshgrid([1:10:190],[1:1:19]);
Ech_elect=l.*(Ee.*(n/k)+Ebf*((n/k)-1)+es*d^2);
Enon_ch_elect=l.*(Ee*(1+k)+(k*Ebf)+es.*((M.*M)/(2*3.14159*k)));
Ech_frame=l.*(((n/k)-m)*Ebf+(el*d^4)+((n/k)-m+1)*Ee);
Enon_ch_frame=l*(Ee+es.*((M.^2)/(2*3.14159*k)));
f1=1./(k.*((n/k)-m+1)); f2=((n/k)-m)./(k.*((n/k)-m+1));
Ech_data=f1.*Nf.*Ech_frame;
Enon_ch_data=f2.*Nf.*Enon_ch_frame;
Ech_iter=Ech_elect+Ech_data;
Enon_ch_iter=Enon_ch_elect+Enon_ch_data;
Estart=(Ech_iter./m)+((n./(k.*m))-1).*(Enon_ch_iter./((n/k)-m));
surf(M,m,Estart);
xlabel('Network Diameter (m)', 'FontSize',14)
ylabel('Head-Sest Size', 'FontSize',14)
zlabel('Estart (J)','FontSize',14)
%>> colormap cool