
%function [acc_temp,acc_ord,acc_inter,acc_tempMean]=tempord1_classify(X,Y,YtempMean,var)
function [zAcc_ovrall, zAcc_ord, zAcc_temp, zAcc_int]=tempord1_classify_RY(X,Y,YtempMean,YcombiCorrected4Temp,var)

nr=6; % runs/nr of scans
Y=Y'; % transpose Y to have a p*n matrix instead of n*p (where n is trials and p is voxels)
YtempMean=YtempMean';

% make vector c with class labels corresponding to Y (72 columns).
append=[1,1,1,2,2,2,3,3,3];
c_temp=repmat(append,1,nr);

append=[1,2,3,1,2,3,1,2,3];
c_ord=repmat(append,1,nr);

append=[1,2,3,4,5,6,7,8,9];
c_ovr=repmat(append,1,nr);

%%%%%Mean T1, T2, T3 averaged over O1, O2, O3 Additional classification analysis
append=[1,2,3];
c_tempMean=repmat(append,1,nr);


%generate (temporarily) a row in Y that indicates which run a trial
%belongs to
run=kron([1:nr],ones(1,9));

%% Generate column indices for cross-validation, where cell i contains
% column indices of the respective test and train set

%%% Classification across all combinations, then correct row and column?
for i=1:nr
    test{i}=find(run==i); % fprintf('test:'); display(test{i});
    train{i}=find(run~=i); % fprintf('train:'); display(train{i});
end;
data=tempord3_combinedclass(Y,c_ovr,run,train,test);
acc_ovr=data(1,1);
% acc_temp=data(2,1);
% acc_ord=data(3,1);

%%Decomposition method
% acc_inter=data(1,:);
% acc_temp=data(2,:);
% acc_ord=data(3,:);

% %%% Classification with a "dedicated" classifier:
% data=tempord3_combinedclass(Y,c_temp,run,train,test);
% acc_tempClass=data(1,1);
% 
% data=tempord3_combinedclass(Y,c_ord,run,train,test);
% acc_ordClass=data(1,1);

%%% Classification across all combinations corrected for Temp patterns
data=tempord3_combinedclass_corrected4Main(Y,c_ovr,run,train,test);
acc_inter=data(1,1);




%%
%%% Classification across mean T1, T2, T3 averaged over O1, O2, O3 Additional classification analysis
run=kron([1:nr],ones(1,3));

for i=1:nr
    test{i}=find(run==i); % fprintf('test:'); display(test{i});
    train{i}=find(run~=i); % fprintf('train:'); display(train{i});
end;
data=tempord3_eventclass(YtempMean,c_tempMean,run,train,test);
acc_tempMean=data(1,1);

%%
%%% Temp One-out classification, e.g. T1, T2, T3 by O1, O2, Tested on T1, T2,
%%% T3 by O3
run=kron([1:nr],ones(1,9));
run=run';
oneout=[1 0 0;0 1 0;0 0 1]; oneout=repmat(oneout,1,3*nr)';
j=0;
for i=1:3:nr*3
    j=j+1;
    test{i}   =find(run==j & oneout(:,1)==1); % fprintf('test:'); display(test{i});
    test{i+1} =find(run==j & oneout(:,2)==1); % fprintf('test:'); display(test{i});
    test{i+2} =find(run==j & oneout(:,3)==1); % fprintf('test:'); display(test{i});
    
    train{i}  =find(run~=j & oneout(:,1)~=1); % fprintf('train:'); display(train{i});
    train{i+1}=find(run~=j & oneout(:,2)~=1); % fprintf('train:'); display(train{i});
    train{i+2}=find(run~=j & oneout(:,3)~=1); % fprintf('train:'); display(train{i});
end;
data=tempord3_eventclass(Y,c_temp,run,train,test);
acc_tempOneOut=data(1,1);

%%%Ord One-out classification
%%% T1_O1 T1_O2 T1_O3 %%% T2_O1 T2_O2 T2_O3 %%% T3_O1 T3_O2 T3_O3
run=kron([1:nr],ones(1,9));
run=run';
oneout=[0 0 1; 0 0 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0];
oneout=repmat(oneout,6,1);
j=0;
for i=1:3:nr*3
    j=j+1;
    test{i}   =find(run==j & oneout(:,1)==1); % Classify Ord with T3
    test{i+1} =find(run==j & oneout(:,2)==1); % Classify Ord with T2
    test{i+2} =find(run==j & oneout(:,3)==1); % Classify Ord with T1
    
    train{i}  =find(run~=j & oneout(:,1)~=1); % Train on Ord with T1 and T2
    train{i+1}=find(run~=j & oneout(:,2)~=1); % Train on Ord with T1 and T3
    train{i+2}=find(run~=j & oneout(:,3)~=1); % Train on Ord with T2 and T3
end;
data=tempord3_eventclass(Y,c_ord,run,train,test);
acc_ordOneOut=data(1,1);

%% Z-transform



%*******For dimensionality analysis
% betas=9;
% categories=9;
% mu=1/categories;
% N=nr*betas;
% sigma=sqrt(mu*(1-mu)*1/N);
% 
% mu=repmat(mu,1,length(acc_inter));%
% sigma=repmat(sigma,1,length(acc_inter));
% zacc_inter=(acc_inter-mu)./sigma;
%*******

%inter
betas=9;
categories=9;
mu=1/categories; %
N=nr*betas;
sigma=sqrt(mu*(1-mu)*1/N);
zAcc_ovrall=(acc_ovr-mu)/sigma;



zAcc_int=(acc_inter-mu)/sigma;


% temp/ord
% betas=9;
% categories=3;
% mu=1/categories; 
% N=nr*betas;
% sigma=sqrt(mu*(1-mu)*1/N);
% zacc_temp=(acc_temp-mu)/sigma;
% zacc_ord=(acc_ord-mu)/sigma;
% 
% zacc_tempClass=(acc_tempClass-mu)/sigma;
% zacc_ordClass=(acc_ordClass-mu)/sigma;



% % tempMean
% betas=3;
% categories=3;
% mu=1/categories;
% N=nr*betas;
% sigma=sqrt(mu*(1-mu)*1/N);
% zacc_tempMean=(acc_tempMean-mu)/sigma;

% timing
betas=3*3; %One-out for training; repeated 3x for each test run
categories=3;
mu=1/categories; 
N=nr*betas;
sigma=sqrt(mu*(1-mu)*1/N);
zAcc_temp=(acc_tempOneOut-mu)/sigma;


% order
betas=3*3; %One-out for training; repeated 3x for each test run
categories=3;
mu=1/categories; 
N=nr*betas;
sigma=sqrt(mu*(1-mu)*1/N);
zAcc_ord=(acc_ordOneOut-mu)/sigma;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % accuracy
% acc_temp=crossval_takemultipleout(@classify_lda_KclassesQuicker,Y,c_temp,train,test);
% acc_ord=crossval_takemultipleout(@classify_lda_KclassesQuicker,Y,c_ord,train,test);
% acc_inter=crossval_takemultipleout(@classify_lda_KclassesQuicker,Y,c_inter,train,test);
% 
% 
% %% Mean T1, T2, T3 averaged over O1, O2, O3 Additional classification analysis
% run=kron([1:nr],ones(1,3));
% 
% for i=1:nr
%     test{i}=find(run==i); % fprintf('test:'); display(test{i});
%     train{i}=find(run~=i); % fprintf('train:'); display(train{i});
% end;
% acc_tempMean=crossval_takemultipleout(@classify_lda_KclassesQuicker,YtempMean,c_tempMean,train,test);
% 
% %disp(acc_temp);
% 
% % z transform
% examples=3;
% mu=0.33;
% N=nr*examples;
% sigma=sqrt(mu*(1-mu)*1/N);
% 
