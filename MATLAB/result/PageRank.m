function PageRank( L, sigma )  
  
%��ȡͼ���ܽ����  
N= size(L,1);   % N���ܽ����(��ҳ����)  
q= 0.85;          % Ĭ������ϵ��  
% ����D  
d= sum(L,2);      % dΪ��ʾ������ȵ�������  
D= diag(d);         
% ����M= L'*D^(-1)   
M= L'*inv(D);      
e= ones(N,1);  
a= (d==0);  % aΪ����������������ҳ����������  
            % ���i��������ȡֵ�ɵ�i����վ�Ƿ�Ϊ��������վ������, ������Ϊ1, ����Ϊ0  
% ����S  
S= M + e*a'/N;     
% �������յĸ���ת�ƾ���  
G= q*S + (1-q)*e*e'/N;       
  
  
  
H0= zeros(N,1);  % ��ʼ��H0  
H1= ones(N,1);   % Ĭ�ϳ�ʼȨ��������Ϊ1  
count= 0;        % ��ʼ����������  
% ������߿�ʼ  
  
while  norm(H1-H0) >= sigma  
    H0= H1;  
    H1= G*H0;  
    count= count+1;  
end  % while������  
  
% �����H1���ս�������,sort()�в���'descend'Ϊ����,'ascend'Ϊ����  
[Rank,index]= sort(H1,'descend');         
% Rank�����б������ս������е�����  
%index�����д洢�Ž���������Ӧ��ԭʼ�����е����  
  
% ��ӡ���  
fprintf('After %d steps PageRank finally converges\n', count)  
fprintf('Rank               Weight              PageNumber\n')  
fprintf('-------------------------------------------------\n')  
for  k= 1:N  
    fprintf('%-3d                %-.7f            %-3d\n',...  
            k, Rank(k), index(k))  
end   % forѭ������  
end    
% PageRank��������  