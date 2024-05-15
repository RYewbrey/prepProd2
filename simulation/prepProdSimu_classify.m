function [zAcc_ord, zAcc_temp, zAcc_int, zAcc_int_subtract]=prepProdSimu_classify(Y)
%outputs are z-scored accuracy

%%%Regular Z-scoring (if 2x2 subtraction issue is ignored in integrated)
numTests=6;
numCat=4;
mu=1/numCat; %mu=0.25;
N=numTests*numCat;
sigma=sqrt(mu*(1-mu)*1/N);

nr=6; % runs/nr of scans
Y=Y'; % transpose Y to have a p*n matrix instead of n*p (where n is trials and p is voxels)

%generate a row that indicates which run a trial
%belongs to
run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9));

% make vector c with class labels corresponding to Y.
c_ovrall=repmat([1:4],1,nr);%ignores unique 2x2 design pattern
c_ovrall_subtract=repmat([1 2 2 1],1,nr);%due to unique subtraction pattern that occurs with a 2x2 design

%%% Generate column indices for cross-validation, where cell i contains
%%% column indices of the respective test and train set
for i=1:nr
    test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
    train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
end

%%% Integrated decoding %%%
acc_int = prepProd2_combinedclass_corrected4Main(Y,c_ovrall,run,train,test);
zAcc_int = (acc_int-mu)/sigma; %z_accuracy=(accuracy-mu)/sigma;

%%%Subtraction Z-scoring (if 2x2 subtraction issue addressed)
numTests=6;
numCat=2;
mu=1/numCat; %mu=0.25;
N=numTests*numCat;
sigmaSubtract=sqrt(mu*(1-mu)*1/N);

%%% Integrated decoding considering subtraction %%%
acc_int_subtract = prepProd2_combinedclass_corrected4Main(Y,c_ovrall_subtract,run,train,test);
zAcc_int_subtract = (acc_int_subtract-mu)/sigmaSubtract; %z_accuracy=(accuracy-mu)/sigma;

%%%One-out Z-scoring
takeOneOutIter=2;
numTests=6;
numCat=2;
mu=1/numCat; %mu=0.5;
N=numTests*numCat*takeOneOutIter;
sigmaOneOut=sqrt(mu*(1-mu)*1/N);

%%%%%%% ORDER DECODING %%%%%%%
c_ord=repmat([1 1 2 2],1,nr); % extract conditions

% Generate column indices for cross-validation, where
% cell i contains column indices of the respective test and
% train set
%%%
oneout=[1 0; 0 1];
oneout=repmat(oneout,1,12);
j = 0;
%trainOrd = cell(max(nr),1); testOrd = cell(max(nr),1);

for i=1:2:nr*2
    j=j+1;
    testOrd{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
    testOrd{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
    
    
    trainOrd{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
    trainOrd{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
    
end

acc_ord = combinedclass(Y,c_ord,run,trainOrd,testOrd);
zAcc_ord = (acc_ord-mu)/sigmaOneOut; % z_accuracy=(accuracy-mu)/sigma;

%%%%%%% TIMING DECODING %%%%%%%
c_temp=repmat([1 2 1 2],1,nr); % extract conditions

% Generate column indices for cross-validation, where
% cell i contains column indices of the respective test and
% train set
%%%
oneout=[1 1 0 0; 0 0 1 1];
oneout=repmat(oneout,1,6);
j=0;
for i=1:2:nr*2
    j=j+1;
    testTemp{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
    testTemp{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
    
    
    trainTemp{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
    trainTemp{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
    
end

acc_temp = combinedclass(Y,c_temp,run,trainTemp,testTemp);
zAcc_temp = (acc_temp-mu)/sigmaOneOut; % z_accuracy=(accuracy-mu)/sigma;


