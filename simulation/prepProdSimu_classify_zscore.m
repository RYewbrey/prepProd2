function [zAcc_ovrall, zAcc_ord, zAcc_temp, zAcc_int, acc_int]=prepProdSimu_classify_zscore(Y)
%Converts all to respective Z scores
%Bugged for integrated - 2x2 design constraints mean that two patterns
%always end the same. Therefore, its chance value is actually 0.50.

nr=6; % runs/nr of scans
Y=Y'; % transpose Y to have a p*n matrix instead of n*p (where n is trials and p is voxels)

numTests=6; %runs
numCat=4; %sequences (within phase)
mu=1/numCat; %mu=0.25;
N=numTests*numCat;
sigma=sqrt(mu*(1-mu)*1/N);

%%%%%%% OVERALL DECODING %%%%%%%
% make vector c with class labels corresponding to Y.
append=[1,2,3,4];
c_ovrall=repmat(append,1,nr);

%generate (temporarily) a row in Y that indicates which run a trial
%belongs to
run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN

%%% Generate column indices for cross-validation, where cell i contains
%%% column indices of the respective test and train set

%%% Classification across all combinations, then correct row and column?
train = cell(1, max(nr)); test = cell(1, max(nr));

for i=1:nr
    test{i}=find(run==i); % fprintf('test:'); display(test{i});
    train{i}=find(run~=i); % fprintf('train:'); display(train{i});
end

%%% Overall decoding %%%
acc_ovrall = combinedclass(Y,c_ovrall,run,train,test);
zAcc_ovrall = (acc_ovrall-mu)/sigma; %z_accuracy=(accuracy-mu)/sigma;

%%% Integrated decoding %%%
acc_int = prepProd2_combinedclass_corrected4Main(Y,c_ovrall,run,train,test);
zAcc_int = (acc_int-mu)/sigma; %z_accuracy=(accuracy-mu)/sigma;

%%%Change Z-scoring to one-out
takeOneOutIter=2;
numTests=6;
numCat=2;
mu=1/numCat; %mu=0.5;
N=numTests*numCat*takeOneOutIter;
sigma=sqrt(mu*(1-mu)*1/N);

%%%%%%% ORDER DECODING %%%%%%%
c_ord=repmat([1 1 2 2],1,nr); % extract conditions

% Generate column indices for cross-validation, where
% cell i contains column indices of the respective test and
% train set
%%%
oneout=[1 0; 0 1];
oneout=repmat(oneout,1,12);
j = 0;
trainOrd = cell(max(nr),1); testOrd = cell(max(nr),1);

for i=1:2:nr*2
    j=j+1;
    testOrd{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
    testOrd{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
    
    
    trainOrd{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
    trainOrd{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
    
end

acc_ord = combinedclass(Y,c_ord,run,trainOrd,testOrd);
zAcc_ord = (acc_ord-mu)/sigma; % z_accuracy=(accuracy-mu)/sigma;

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
zAcc_temp = (acc_temp-mu)/sigma; % z_accuracy=(accuracy-mu)/sigma;



