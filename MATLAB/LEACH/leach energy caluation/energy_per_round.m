el=0.0013*10^-12; es=10*10^-12; Ee=50*10^-9; Ebf=5*10^-9;
n=1000; k=50; d=100; l=4000; Nf=10000; m=1;
[M,k]=meshgrid([1:20:200],[5:5:50]);
Ech_elect=l.*(Ee.*(n./k)+Ebf.*((n./k)-1)+es*d^2);
Enon_ch_elect=l.*(Ee.*(1+k)+(k*Ebf)+es.*((M.*M)./((2*3.14159).*k)));
Enon_ch_frame=l.*(Ee+es.*((M.^2)./((2*3.14159).*k)));
Ech_frame=l.*(((n./k)-m)*Ebf+(el*d^4)+((n./k)-m+1)*Ee);
f1=1./(k.*((n./k)-m+1)); f2=((n./k)-m)./(k.*((n./k)-m+1));
Ech_data=f1.*Nf.*Ech_frame;
Enon_ch_data=f2.*Nf.*Enon_ch_frame;
Ech_iter=Ech_elect+Ech_data;
Enon_ch_iter=Enon_ch_elect+Enon_ch_data;
Estart=(Ech_iter./m)+((n./(k.*m))-1).*(Enon_ch_iter./((n./k)-m));
surf(M,k,Estart);
hold on
m=3;
Ech_elect=l.*(Ee.*(n./k)+Ebf.*((n./k)-1)+es*d^2);
Enon_ch_elect=l.*(Ee.*(1+k)+(k*Ebf)+es.*((M.*M)./((2*3.14159).*k)));
Enon_ch_frame=l.*(Ee+es.*((M.^2)./((2*3.14159).*k)));
Ech_frame=l.*(((n./k)-m)*Ebf+(el*d^4)+((n./k)-m+1)*Ee);
f1=1./(k.*((n./k)-m+1)); f2=((n./k)-m)./(k.*((n./k)-m+1));
Ech_data=f1.*Nf.*Ech_frame;
Enon_ch_data=f2.*Nf.*Enon_ch_frame;
Ech_iter=Ech_elect+Ech_data;
Enon_ch_iter=Enon_ch_elect+Enon_ch_data;
Estart=(Ech_iter./m)+((n./(k.*m))-1).*(Enon_ch_iter./((n./k)-m));
surf(M,k,Estart);
xlabel('Network Diameter (m)', 'FontSize',14)
ylabel('# Of Clusters', 'FontSize',14)
zlabel('Estart (J)','FontSize',14)
%>> colormap cool