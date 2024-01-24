function [C,aver_C]=Degree_Clustering_Coefficient(A)
%% ������ͼ�и��ڵ�ľ���ϵ������������ľ���ϵ��
%% ����㷨�����ÿ���ڵ�ľ���ϵ������ĳ�ڵ�������ھӣ���Щ�ھӽڵ㹹��һ����ͼ
%% ��A�г������ͼ���ڽӾ��󣬼�����ͼ�ı������ٸ��ݾ���ϵ���Ķ��壬��������ýڵ�ľ���ϵ��
%A��������������������ͼ���ڽӾ���
%C��������������������ͼ���ڵ�ľ���ϵ��
%aver����������������������ͼ�ľ���ϵ��
N=size(A,2);
C=zeros(N,1);
for i=1:N
    aa=find(A(i,:)==1); %Ѱ����ͼ���ھӽڵ�
    if isempty(aa)
    disp(['�ڵ�',int2str(i),'Ϊ�����ڵ㣬�����ϵ����ֵΪ0']);
    C(i)=0;
    else
    m=length(aa); % Calculate the degree of the graph
    if m==1
    disp(['�ڵ�',int2str(i),'ֻ��һ���ھӽڵ㣬�����ϵ����ֵΪ0']);
    C(i)=0;
    else
    B=A(aa,aa); % ��ȡ��ͼ���ڽӾ���
    C(i) = 2*length(find(B==1))/(m*(m-1));
    % Calculate the importance of the node
    end
    end
end
aver_C=mean(C);
    