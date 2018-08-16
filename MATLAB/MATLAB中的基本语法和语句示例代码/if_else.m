function if_else(a)

if mod(a,2)==0
    disp('a为2的倍数');
elseif mod(a,3)==0
    disp('a为3的倍数');
elseif mod(a,5)==0
    disp('a为5的倍数');
else
    disp('a不是2、3、5的倍数');
end
