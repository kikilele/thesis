function son_function( )        %����������λ�����Ϸ�
%�Ӻ�������
max1=find_max(1,2,3)
max2=find_max(7,3,9)
max3=find_max(2,10,2)
max4=find_max(1,2,2)
max5=find_max(4,3,4)


function max=find_max(a,b,c)    %�Ӻ���
if (a>=b)&&(a>=c)
    max=a;
elseif (b>=a)&&(b>=c)
    max=b;
else
    max=c;
end
