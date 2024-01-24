clear all
clc

addpath('.\LMCI\LMCIAAL\3ROISignals_FunImgARWSDCF');

%% Read data
LMCINodeCC3 = zeros(90,31);
for num = 28:31
    if num > 9
        load(['ROICorrelation_FisherZ_Sub_0',num2str(num),'.mat'])
    else 
        load(['ROICorrelation_FisherZ_Sub_00',num2str(num),'.mat'])
    end
    %% ȥ��С�Բ������ڽӾ���
    A90=zeros(90,90); % Ȩ�ؾ���
    H=zeros(90,90); % 0-1 �ڽӾ���
    for i = 1:90
          for j = 1:90
              if ROICorrelation_FisherZ(i,j) < 0 || i==j
                 A90(i,j) = 0;
              else 
                 A90(i,j) = ROICorrelation_FisherZ(i,j);
              end
              if A90(i,j) > 0
                  H(i,j) = 1;
              else 
                  H(i,j) = 0;
              end
          end
    end
    %% ÿ��������Ϊһ���ڵ㣬��ȡ�ڵ�����
    % Clustering coefficient
    [C,aver_C]=Degree_Clustering_Coefficient(H);
    % pagerank algorithm
    P = zeros(90,1);
    [m,n] = size(H);
    for i = 1:n
        sa(i) = sum(A90(i,:),2);
        Na(i) = sum(H(i,:),2); %���������A���еĺ�
        %d = 0.85; %��������
        P(i) = sa(i)/Na(i);
        
        finalIndex(i) = 0.5*C(i) + 0.5*P(i);
 
    end
     
    LMCINodeCC3(:,num) = finalIndex;
    save('LMCINodeCC3.mat');
end
    