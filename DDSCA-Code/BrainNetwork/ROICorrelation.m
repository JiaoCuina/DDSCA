% addpath('.\ROICorrelation');
addpath('.\LMCI\LMCIbrain\3ROISignals_FunImgARWSDCF');

%% 去除小脑，90*90
for num = 28:31
    if num > 9
        load(['ROICorrelation_FisherZ_Sub_0',num2str(num),'.mat'])
    else 
        load(['ROICorrelation_FisherZ_Sub_00',num2str(num),'.mat'])
    end

    %% 提取上对角元素，转化成列向量
        CorrelationDiag=triu(ROICorrelation_FisherZ)-diag(diag(ROICorrelation_FisherZ));

        length=size(CorrelationDiag,1);
        %  控制A列变量
        m = 2;
        %  控制B列变量
        n = 1;
        for i =1:length
            for j =m:length
                LMCICorrelationVector3{num}(n,1)=CorrelationDiag(i,j); 
                n=n+1;
            end
            m=m+1;
        end
     save('LMCICorrelationVector3.mat');
end

