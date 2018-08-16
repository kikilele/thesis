clear
clc
%�������ɿռ�TDOA���棨ţ�ٵ�������
%ע������Ϊ�����Ƚ��յ��źŵĻ�վ��Ϊ0������Ϊ����ԭ�㽨������ϵ

%������ʼ��
star_x=[0 0 20 20 ];%��λkm
star_y=[0 20 20 0 ];
star_z=[-0.3 0 0 0 ];

xyz=[200 0 10];%��λkm

d0=sqrt((star_x(1)-xyz(1)).^2+(star_y(1)-xyz(2)).^2+(star_z(1)-xyz(3)).^2);
d1=sqrt((star_x(2)-xyz(1)).^2+(star_y(2)-xyz(2)).^2+(star_z(2)-xyz(3)).^2);
d2=sqrt((star_x(3)-xyz(1)).^2+(star_y(3)-xyz(2)).^2+(star_z(3)-xyz(3)).^2);
d3=sqrt((star_x(4)-xyz(1)).^2+(star_y(4)-xyz(2)).^2+(star_z(4)-xyz(3)).^2);

c=3e8;
Tmiss=20*c*1e-9*1e-3;%��λns�������Ǿ��뵥λΪkm
delta_r1=d1-d0+Tmiss*rand(1,1);
delta_r2=d2-d0+Tmiss*rand(1,1);
delta_r3=d3-d0+Tmiss*rand(1,1);

%ģ�ͷ��̽ⷨ����һ�ֿռ�ʱ�λ�����㷨��
syms x y z
r0=sqrt((x-star_x(1)).^2+(y-star_y(1)).^2+(z-star_z(1)).^2);
r1=sqrt((x-star_x(2)).^2+(y-star_y(2)).^2+(z-star_z(2)).^2);
r2=sqrt((x-star_x(3)).^2+(y-star_y(3)).^2+(z-star_z(3)).^2);
r3=sqrt((x-star_x(4)).^2+(y-star_y(4)).^2+(z-star_z(4)).^2);

F=[r1-r0-delta_r1;
   r2-r0-delta_r2;
   r3-r0-delta_r3];

F_dif=[(x-star_x(2))./r1-(x-star_x(1))./r0, (y-star_y(2))./r1-(y-star_y(1))./r0, (z-star_z(2))./r1-(z-star_z(1))./r0;
       (x-star_x(3))./r2-(x-star_x(1))./r0, (y-star_y(3))./r2-(y-star_y(1))./r0, (z-star_z(3))./r2-(z-star_z(1))./r0;
       (x-star_x(4))./r3-(x-star_x(1))./r0, (y-star_y(4))./r3-(y-star_y(1))./r0, (z-star_z(4))./r3-(z-star_z(1))./r0];

F_newton=inv(F_dif)*F;

start=[220;5;12];
result=Newton(start,F_newton,0.3)

GDOP=(xyz(1)-result(1)).^2+(xyz(2)-result(2)).^2+(xyz(3)-result(3)).^2;
disp([sqrt(GDOP)]);
