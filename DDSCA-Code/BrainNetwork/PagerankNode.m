addpath('.\LMCI\LMCIbrain\3ROISignals_FunImgARWSDCF');

%% Read data
LMCINode3 = zeros(246,31);
for num = 28:31
    if num > 9
        load(['ROICorrelation_FisherZ_Sub_0',num2str(num),'.mat'])
    else 
        load(['ROICorrelation_FisherZ_Sub_00',num2str(num),'.mat'])
    end
    %% 去除小脑
    A246=zeros(246,246); % 权重矩阵
    H=zeros(246,246); % 0-1 邻接矩阵
    for i = 1:246
          for j = 1:246
              if ROICorrelation_FisherZ(i,j) < 0 || i==j
                 A246(i,j) = 0;
              else 
                 A246(i,j) = ROICorrelation_FisherZ(i,j);
              end
              if A246(i,j) > 0
                  H(i,j) = 1;
              else 
                  H(i,j) = 0;
              end
          end
    end
    
    P = zeros(246,1);
    [m,n] = size(H);
    for i = 1:n
        sa(i) = sum(A246(i,:),2);
        Na(i) = sum(H(i,:),2); %按列求矩阵A各行的和
        %d = 0.85; %阻尼因子

        P(i) = sa(i)/Na(i);
    end
    LMCINode3(:,num) = P;
    save('LMCINode3.mat');
end