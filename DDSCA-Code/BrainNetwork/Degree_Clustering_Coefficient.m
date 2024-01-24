function [C,aver_C]=Degree_Clustering_Coefficient(A)
%% 求网络图中各节点的聚类系数及整个网络的聚类系数
%% 求解算法：求解每个节点的聚类系数，找某节点的所有邻居，这些邻居节点构成一个子图
%% 从A中抽出该子图的邻接矩阵，计算子图的边数，再根据聚类系数的定义，即可算出该节点的聚类系数
%A————————网络图的邻接矩阵
%C————————网络图各节点的聚类系数
%aver———————整个网络图的聚类系数
N=size(A,2);
C=zeros(N,1);
for i=1:N
    aa=find(A(i,:)==1); %寻找子图的邻居节点
    if isempty(aa)
    disp(['节点',int2str(i),'为孤立节点，其聚类系数赋值为0']);
    C(i)=0;
    else
    m=length(aa); % Calculate the degree of the graph
    if m==1
    disp(['节点',int2str(i),'只有一个邻居节点，其聚类系数赋值为0']);
    C(i)=0;
    else
    B=A(aa,aa); % 抽取子图的邻接矩阵
    C(i) = 2*length(find(B==1))/(m*(m-1));
    % Calculate the importance of the node
    end
    end
end
aver_C=mean(C);
    