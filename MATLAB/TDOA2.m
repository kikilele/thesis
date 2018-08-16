clear
clc
%�������ɿռ�TDOA���棨��������
%ע������Ϊ�����Ƚ��յ��źŵĻ�վ��Ϊ0������Ϊ����ԭ�㽨������ϵ

%������ʼ��
star_x=[0e-3    0e-3  1000e-3 1000e-3 ];    %��λkm
star_y=[0e-3    1000e-3  1000e-3 0e-3 ];
star_z=[-300e-3 0e-3  0e-3  0e-3 ];
x=10*(-20:20);
y=10*(-20:20);
z=10*ones(1,length(x));
xyz_all=[x ;y ;z];%��λkm
GDOP_avr=zeros(length(x),length(y));
clear x y z
%----------------------------��1��2��ѭ����������ά����
for row=1:size(xyz_all,2)
for column=1:size(xyz_all,2)

xyz=[xyz_all(1,row),xyz_all(2,column),xyz_all(3,1)];

d0=sqrt((star_x(1)-xyz(1)).^2+(star_y(1)-xyz(2)).^2+(star_z(1)-xyz(3)).^2);
d1=sqrt((star_x(2)-xyz(1)).^2+(star_y(2)-xyz(2)).^2+(star_z(2)-xyz(3)).^2);
d2=sqrt((star_x(3)-xyz(1)).^2+(star_y(3)-xyz(2)).^2+(star_z(3)-xyz(3)).^2);
d3=sqrt((star_x(4)-xyz(1)).^2+(star_y(4)-xyz(2)).^2+(star_z(4)-xyz(3)).^2);

c=3e8;
Tmiss=5*c*1e-9*1e-3;%��λns�������Ǿ��뵥λΪkm
%----------------------------��3��ѭ��������ĳ����ľ���
N=20;
count=N;%�������ɹ��Ĵ������Ա�ȡƽ��
x_N=0;y_N=0;z_N=0;
for i=1:N
delta_r1=d1-d0+Tmiss*rand(1,1);
delta_r2=d2-d0+Tmiss*rand(1,1);
delta_r3=d3-d0+Tmiss*rand(1,1);

%ģ�ͷ��̽ⷨ������վ��Դʱ�λ������
A=[star_x(1)-star_x(2),star_y(1)-star_y(2),star_z(1)-star_z(2);
   star_x(1)-star_x(3),star_y(1)-star_y(3),star_z(1)-star_z(3);
   star_x(1)-star_x(4),star_y(1)-star_y(4),star_z(1)-star_z(4)];

k=0.5*[delta_r1.^2+star_x(1).^2+star_y(1).^2+star_z(1).^2-star_x(2).^2-star_y(2).^2-star_z(2).^2;
       delta_r2.^2+star_x(1).^2+star_y(1).^2+star_z(1).^2-star_x(3).^2-star_y(3).^2-star_z(3).^2;
       delta_r3.^2+star_x(1).^2+star_y(1).^2+star_z(1).^2-star_x(4).^2-star_y(4).^2-star_z(4).^2];
   
A=inv(A);
   
m=[A(1,1)*k(1)+A(1,2)*k(2)+A(1,3)*k(3);
   A(2,1)*k(1)+A(2,2)*k(2)+A(2,3)*k(3);
   A(3,1)*k(1)+A(3,2)*k(2)+A(3,3)*k(3)];

n=[A(1,1)*delta_r1+A(1,2)*delta_r2+A(1,3)*delta_r3;
   A(2,1)*delta_r1+A(2,2)*delta_r2+A(2,3)*delta_r3;
   A(3,1)*delta_r1+A(3,2)*delta_r2+A(3,3)*delta_r3];

a=n(1)^2+n(2)^2+n(3)^2-1;
b=(m(1)-star_x(1))*n(1)+(m(2)-star_y(1))*n(2)+(m(3)-star_z(1))*n(3);
c=(m(1)-star_x(1))^2+(m(2)-star_y(1))^2+(m(3)-star_z(1))^2;

if b^2-a*c < 0
    count=count-1;
else
    r0=[(-b+sqrt(b^2-a*c))/a,(-b-sqrt(b^2-a*c))/a];
    x=m(1)+n(1)*r0;
    y=m(2)+n(2)*r0;
    z=m(3)+n(3)*r0;
end 

x_N=x+x_N;
y_N=y+y_N;
z_N=z+z_N;
% GDOP_avr(row,column)=GDOP_avr(row,column)+GDOP;
end

    x=x_N/count;
    y=y_N/count;
    z=z_N/count;
    GDOP1=(xyz(1)-x(1)).^2+(xyz(2)-y(1)).^2+(xyz(3)-z(1)).^2;
    GDOP2=(xyz(1)-x(2)).^2+(xyz(2)-y(2)).^2+(xyz(3)-z(2)).^2;
    GDOP=sqrt(min(GDOP1,GDOP2));%ȡ��С�ļ����Ѿ����ģ������
    if count==0;
        GDOP=25;
    end
GDOP_avr(row,column)=GDOP;
end
end

meshc(10*(-20:20),10*(-20:20),100*GDOP_avr);
zmax=max(max(GDOP_avr));
zmin=min(min(GDOP_avr));
n0=(zmax-zmin)/45;
for nn0=zmin:n0:zmax
  [cs,h]=contour(10*(-20:20),10*(-20:20),GDOP_avr,nn0);
  hold on
end