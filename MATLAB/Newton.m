%牛顿迭代法解非线性方程组
%输入参数：迭代初值，迭代方程，精度要求
function [h]=Newton(start_newton,F_newton,p_newton)
F_temp=start_newton;
temp=[0;0;0];   

while sum(abs(temp-F_temp))>p_newton
temp=F_temp;    
F_temp=F_temp-subs(F_newton,{'x','y','z'},{F_temp(1),F_temp(2),F_temp(3)});
end

h=F_temp;