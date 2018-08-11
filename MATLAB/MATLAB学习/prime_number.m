clear
sum=5;         %求0～100素数之和
ss=0;          %用来标定是否是素数，0表示不是
prime=[2 3];     %用来存放素数，2，3为素数，先放置在prime矩阵中
for i=4:100
    for j=2:fix(sqrt(i))
        if mod(i,j)==0
            ss=0;     %能被整除，说明i不是素数，用ss=0来表示
            break;    %能被整除，跳出内循环
        else 
            ss=1;
        end
    end
    if ss==1          %是素数，保存至prime矩阵，并求和
        prime=[prime,i];
        sum=sum+i;
    end
end
fprintf('%d,',prime);
fprintf('\nThe prime''s sum is %d',sum);

        