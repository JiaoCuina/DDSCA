clear all
clc

%% importing data and making network parameters in an optimum state and partition Data for the 5 fold test.
load Brain246.mat
% load Brain90.mat
% load img_name2.mat
% load snp_name2.mat
%% Use the neural network parameters below.
%randseed=8409;
rcov1=1e-4; rcov2=1e-4;
% Hidden activation type.
hiddentype='sigmoid';
% Architecture (hidden layer sizes) for genotypes data neural network.

NN1=[512 512 512 512 85]; 
%  Architecture (hidden layer sizes)  for brain network phenotypes data neural network.
NN2=[512 512 512 512 30135];
     
% Weight decay parameter.
l2penalty=1e-4;

%% Run Deep Neural Network with SGD.
% Minibatchsize.
batchsize=25;
% Learning rate.
eta0=0.01;
% Rate in which learning rate decays over iterations.
% 1 means constant learning rate.
decay=1;
% Momentum.
momentum=0.99;
% How many passes of the data you run SGD with.
maxepoch=200;
addpath ./deepnet/

X1 = SNPdata2; X2 = BrainNetCC; label = Score2(:,5);
[N,D1]=size(X1); [~,D2]=size(X2);
%% Relief
[ranks,weights] = relieff(X2,label,10);
[W,index] = sort(weights,'descend');
 %% Sort
  Edge = zeros(N,D2);
  for i = 1:size(ranks,2)
      Edge(:,i)=X2(:,ranks(i));
  end
  
  % j = 3000:300:6000 找最有特征值
for j = 4200
      YY = zeros(N,j);
      YY(:,:) = Edge(:,1:j);
      BB(:,:) = [NodeCCbrain, YY];
      [~,D3]=size(BB);
end
X3=BB; 
  %% Set genotypes data architecture.
  Layersizes1=[D1 NN1];  Layertypes1={};
  for nn1=1:length(NN1)-1;
    Layertypes1=[Layertypes1, {hiddentype}];
  end
  % I choose to set the last layer to be linear.
  Layertypes1{end+1}='linear';
  %% Set brain network phenotypes data architecture.
 Layersizes2=[D3 NN2];  Layertypes2={};
 for nn2=1:length(NN2)-1;
   Layertypes2=[Layertypes2, {hiddentype}];
 end
Layertypes2{end+1}='linear';
  %% Random initialization of weights.
  F1=deepnetinit(Layersizes1,Layertypes1);
  F2=deepnetinit(Layersizes2,Layertypes2);
  
  for j=1:length(F1)  F1{j}.l=l2penalty;  end
  for j=1:length(F2)  F2{j}.l=l2penalty;  end
  
  %% the outputs at the top layer.
  FX1=deepnetfwd(X1,F1); FX2=deepnetfwd(X3,F2);

%the self-representation matrix is learned for reconstructing the source data at the top layer.
 
Q = [0.1];
for num = 1:length(Q)
     
            options.lambda = Q(num);
            opts = [];
            opts.init = 0;
            opts.tFlag = 10; opts.maxIter = 100;
            opts.rFlag = 10^-5;
            opts.rsL2 = 0;
            options.opts = opts;
            options.label = Score2(:,5);
            label = Score2(:,5);
% 
           wSNPdata = f_SR(FX1', options);
           wSNPdata=wSNPdata+wSNPdata';
     
        
           wBNdata = f_SR(FX2', options);
           wBNdata=wBNdata+wBNdata';
           
         
            sSNPdata = wSNPdata*X1;
            sBNdata = wBNdata*X3;
%% Kfold cross validation

    Z = label;

    kk=1; 
    kfold = 5;
    [tcv,fcv]=f_myCV(label',kfold,kk);
    for cc = 1:kfold
        trLab=tcv{cc}';
        teLab=fcv{cc}';
        X{1,cc}=sSNPdata(trLab,:);  
        Y{1,cc}=sBNdata(trLab,:);
        Label{1,cc}= Z(trLab,:);

        Xt{1,cc}=sSNPdata(teLab,:);
        Yt{1,cc}=sBNdata(teLab,:);
        Labelt{1,cc}=Z(teLab,:);

    end
    
%%  五折交叉验证
    for n=1:5
    
    M1_te=Xt{1,n};
    M2_te=Yt{1,n};
    M3_te=Labelt{1,n};
    M1_tr=X{1,n};
    M2_tr=Y{1,n};
    M3_tr=Label{1,n};
 
    %%
    test_X=M1_te;
    test_Y=M2_te;
    test_Z=M3_te;
    train_X=M1_tr;
    train_Y=M2_tr; 
    train_Z=M3_tr; 
    gnd=Label{1,n};
    %%main regression function
    % Set parameters
    paraset=[0.001 0.01 0.1 1 10 100 1000];
    if n==1
    for i=1:length(paraset)
        para.lambda1 = paraset(i);
        for j=1:length(paraset)
            para.lambda2 = paraset(j);
            kfold=5;
            kk=2;     
            % construct the index of cross_validation for each task.               
            [tcv, fcv]=f_myCV(gnd',kfold,kk); 
            %% begin to 5-fold.

            for cc=1:kfold 
                trLab=tcv{cc}';
                teLab=fcv{cc}';
                X_tr=train_X(trLab,:);
                Y_tr=train_Y(trLab,:);
                Z_tr=train_Z(trLab,:);
                X_te=train_X(teLab,:);
                Y_te=train_Y(teLab,:);
                Z_te=train_Z(teLab,:);
                [w1opt,w2opt,obj] = f_sCCALR(X_tr,Y_tr,Z_tr,para);
               
                %%
                comp_Xopt = X_te * w1opt;
                comp_Yopt = Y_te * w2opt;
                CCopt(cc) = corr(comp_Xopt, comp_Yopt); 
            end
            res_CC(i,j)=mean(CCopt);
        end
    end
 
    tempCC=0;
    for ii=1:length(paraset)
        for jj=1:length(paraset)
           if  res_CC(ii,jj)>tempCC
             tempCC=res_CC(ii,jj);
             paraOpt=[ii,jj];
           end
        end
    end
    paraOpt
    paraFinal.lambda1=paraset(paraOpt(1));
    paraFinal.lambda2=paraset(paraOpt(2));
    end
    
  %%  Solve the Sample imbalance
  [train_X1, train_Y1, train_Z1] = do_oversample(train_X, train_Y, train_Z); 
  [w1,w2,obj] = f_sCCALR(train_X1,train_Y1,train_Z1,paraFinal);

    %%
    comp_Xtr = train_X * w1;
    comp_Ytr = train_Y * w2;
    CCtr(n) = corr(comp_Xtr, comp_Ytr);
    comp_Xte = test_X * w1;
    comp_Yte = test_Y * w2;
    CCte(n) = corr(comp_Xte, comp_Yte);
    Weight1{n}=abs(w1); %SNP
    Weight2{n}=abs(w2); %Brain Phenotype
     
    end
   
res_CCte(num,1)=mean(CCte);
res_CCtr(num,1)=mean(CCtr);
std_CCte(num,1)=var(CCte);
std_CCtr(num,1)=var(CCtr);

 end
%% results selected top SNP features
% % FC
Weight_tr1 = (Weight1{1}+Weight1{2}+Weight1{3}+Weight1{4}+Weight1{5})/5;
top_K = 10;
Robust_selected1 = get_top_K_features(Weight_tr1, top_K);
% % results selected top Node features
% % Node/FC
Weight_tr2=(Weight2{1}+Weight2{2}+Weight2{3}+Weight2{4}+Weight2{5})/5;
top_K = 10;
Robust_selected2 = get_top_K_features(Weight_tr2, top_K);

%% ---------- Figures -----------
% SNP features
figure(1);
hold on;
subplot(414);
colormap('Jet');
imagesc(Weight_tr1', [-0.2 0.2]);%以图像形式形式矩阵
set(gca, 'XTick',1:85,'XTickLabel',snp_name2,'FontSize',10,'FontName','Times New Roman','XTickLabelRotation',90);
set(gca, 'YTick', 1 : 1, 'YTickLabel', {'AD'});
colorbar;

% Brain features
figure(2);
hold on;
subplot(414);
colormap('Jet');
imagesc(Weight_tr2', [-0.2 0.2]);
set(gca, 'XTick',1:90,'XTickLabel',img_name2,'FontSize',10,'FontName','Times New Roman','XTickLabelRotation',90);
set(gca, 'YTick', 1 : 1, 'YTickLabel', {'AD'});
colorbar;

