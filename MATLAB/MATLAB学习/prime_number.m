clear
sum=5;         %��0��100����֮��
ss=0;          %�����궨�Ƿ���������0��ʾ����
prime=[2 3];     %�������������2��3Ϊ�������ȷ�����prime������
for i=4:100
    for j=2:fix(sqrt(i))
        if mod(i,j)==0
            ss=0;     %�ܱ�������˵��i������������ss=0����ʾ
            break;    %�ܱ�������������ѭ��
        else 
            ss=1;
        end
    end
    if ss==1          %��������������prime���󣬲����
        prime=[prime,i];
        sum=sum+i;
    end
end
fprintf('%d,',prime);
fprintf('\nThe prime''s sum is %d',sum);

        