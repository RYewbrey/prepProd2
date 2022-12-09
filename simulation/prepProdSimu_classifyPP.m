function [acc_ovrall, acc_int, acc_ord, acc_temp]=prepProdSimu_classifyPP(Y)

nr=6; % runs/nr of scans
Y=Y'; % transpose Y to have a p*n matrix instead of n*p (where n is trials and p is voxels)

%%%%%%% OVERALL DECODING %%%%%%%
% make vector c with class labels corresponding to Y.
append=[1,2,3,4];
c_ovrall=repmat(append,1,nr);

%generate (temporarily) a row in Y that indicates which run a trial
%belongs to
run=kron(1:nr,ones(1,4));

%%% Generate column indices for cross-validation, where cell i contains
%%% column indices of the respective test and train set

%%% Classification across all combinations, then correct row and column?
train = cell(max(nr),1); test = cell(max(nr),1);

for i=1:nr
    test{i}=find(run==i); % fprintf('test:'); display(test{i});
    train{i}=find(run~=i); % fprintf('train:'); display(train{i});
end;

%%% Overall decoding %%%
acc_ovrall = combinedclass(Y,c_ovrall,run,train,test);

%%% Integrated decoding %%%
acc_int = prepProd2_combinedclass_corrected4Main(Y,c_ovrall,run,train,test);

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
    
end;

acc_ord = combinedclass(Y,c_ord,run,trainOrd,testOrd);



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
    
end;



acc_temp = combinedclass(Y,c_temp,run,trainTemp,testTemp);



