N=size(adj,2);
 D=adj;
 D(D==0)=inf;    %���ڽӾ����Ϊ�ڽӾ�����������ޱ�����ʱ��ֵΪinf����������ľ���Ϊ0.
 for i=1:N           
     D(i,i)=0;       
 end   
 for k=1:N            %Floyd�㷨��������������̾���
     for i=1:N
         for j=1:N
             if D(i,j)>D(i,k)+D(k,j)
                D(i,j)=D(i,k)+D(k,j);
             end
         end
     end
 end
 CC=zeros(1,N);
 D(D==inf)=0;
 for i=1:N
     CC(i)=(N-1)/sum(D(i,:));%Closeness centrality
     %if sum(D(i,:))==inf
        % CC(i)=0;
 %end
 end
  aver_CC=sum(CC)/N;  %ƽ��cc