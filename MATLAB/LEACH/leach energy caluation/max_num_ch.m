el=0.0013*10^-12; es=10*10^-12; Ee=50*10^-9; Ebf=5*10^-9;
n=1000; M=100; d=100; l=4000; Nf=10000;
[k]=[1:1:20];
m=1;
Ech_elect=l.*(Ee.*((n./k)-m+1)+es*d^2);
Enon_ch_elect=l.*(Ee.*(1+k)+es*d^2);
Ech_frame=l.*(((n./k)-m+1)*Ebf+el*d^4+((n./k)-m+1)*Ee);
Enon_ch_frame=l.*(Ee+es*((M^2)./((2*3.14159).*k)));
Ech_data=(1./((n./k)-m+1)).*(Nf./k).*Ech_frame;
Enon_ch_data=(((n./k)-m)./((n./k)-m+1)).*(Nf./k).*Enon_ch_frame;
Ech_iter=Ech_elect+Ech_data;
Enon_ch_iter=Enon_ch_elect+Enon_ch_data;
Eround=(Ech_iter/m)+(((n./(k*m))-1).*Enon_ch_iter)./((n./k)-m);
plot(k,Eround);
hold on
m=3;
Ech_elect=l.*(Ee.*((n./k)-m+1)+es*d^2);
Enon_ch_elect=l.*(Ee.*(1+k)+es*d^2);
Ech_frame=l.*(((n./k)-m+1)*Ebf+el*d^4+((n./k)-m+1)*Ee);
Enon_ch_frame=l.*(Ee+es*((M^2)./((2*3.14159).*k)));
Ech_data=(1./((n./k)-m+1)).*(Nf./k).*Ech_frame;
Enon_ch_data=(((n./k)-m)./((n./k)-m+1)).*(Nf./k).*Enon_ch_frame;
Ech_iter=Ech_elect+Ech_data;
Enon_ch_iter=Enon_ch_elect+Enon_ch_data;
Eround=(Ech_iter/m)+(((n./(k*m))-1).*Enon_ch_iter)./((n./k)-m);
plot(k,Eround);
xlabel('Number Of Clusters', 'FontSize',14)
ylabel('Energy (J)', 'FontSize',14)