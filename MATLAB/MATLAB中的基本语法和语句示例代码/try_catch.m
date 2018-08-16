clear
a=magic(4); %4*4矩阵，行，列，对角线的数字之和均相等
b=eye(3);%单位阵，3×3
try
    c=a*b;        %执行该语句段出现错误，转而执行catch之后的语句段
catch ME
    disp('内部矩阵维度不一致');
end

