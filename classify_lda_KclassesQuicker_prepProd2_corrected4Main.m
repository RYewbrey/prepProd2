function [cpred,Sw,L,K,muK] = classify_lda_KclassesQuicker_prepProd2_corrected4Main(xtrain, ctrain, xtest,regularization)
% Multi-class classification using linear discriminant analysis 
% without prior!!!!
% INPUT:
%   xtrain : training set, p*ctr matrix with c datapoints in p dimensions
%   ctrain : 1*ctr vector with class labels corresponding to xtrain.
%            Only two different class labels can be used.
%   xtest  : test set, p*cte matrix with c2 datapoints in p dimensions 
%   OPTIONS: 
%       'regularization',0.03 
% OUTPUT:    
%   cpred  : 1*cte vector with predicted class labels for xtest.
%
% JD June 2010
% KK Adapted for main factor correction (combined classifier)
if (nargin<4) 
    regularization=0.01;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. TRAIN: Subtract out the mean of T1, T2 and S1, S2 from xtrain
Y=xtrain';
%Temp
k=0;
n=size(Y,1);
for i=1:2:n; %2
    k=k+1;
    YspatMean(k,:)=mean(Y(i:(i+1),:),1); %1
end;
YspatMean2Subtract=kron(YspatMean,[1;1]); %triple mean of voxels

%Ord %%%%
YtempMean2Subtract=[];
k=0;
j=2; %distance between betas/trials with the same timing (T1,T2)
for i=1:4:n; %for each run
    k=k+1;
    YtempMean(k,:)=mean(Y(i:j:(i+1*j),:),1);
    k=k+1;
    YtempMean(k,:)=mean(Y((i+1):j:((i+1)+1*j),:),1);

    YtempMean2SubtractRun=repmat(YtempMean,j,1);
    YtempMean2Subtract=[YtempMean2Subtract;YtempMean2SubtractRun];
    YtempMean2SubtractRun=[]; YtempMean=[]; k=0;
end;

%%Subtract mean T1, T2, T2, from Y, respectively
YcombiCorrected4Main=Y-(YspatMean2Subtract+YtempMean2Subtract);
xtrain=YcombiCorrected4Main';

% Empty Means:
YtempMean=[];
YspatMean=[];
YspatMean2Subtract=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. TEST: Subtract out the mean of T1, T2, T3 and O1, O2, O3 from xtest
%Temp
YspatMean2Subtract=[];
Y=xtest';
k=0;
n=size(Y,1);
for i=1:2:n;
    k=k+1;
    YspatMean(k,:)=mean(Y(i:(i+1),:),1);
end;
YspatMean2Subtract=kron(YspatMean,[1;1]); %triple mean of voxels

%Ord %%%%
YtempMean2Subtract=[];
k=0;
j=2; %distance between betas/trials with the same order (O1,O2,O3)
for i=1:4:n; %for each run
    k=k+1;
    YtempMean(k,:)=mean(Y(i:j:(i+1*j),:),1);
    k=k+1;
    YtempMean(k,:)=mean(Y((i+1):j:((i+1)+1*j),:),1);

    YtempMean2SubtractRun=repmat(YtempMean,j,1);
    YtempMean2Subtract=[YtempMean2Subtract;YtempMean2SubtractRun];
    YtempMean2SubtractRun=[]; YtempMean=[]; k=0;
end;

%%Subtract mean T1, T2, T2, from Y, respectively
YcombiCorrected4Main=Y-(YspatMean2Subtract+YtempMean2Subtract);
xtest=YcombiCorrected4Main';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[P,N]=size(xtrain);     % size of training set 
classes=1:max(ctrain);  % classses we do classification on
cc=size(classes,2);     % class count

muK=zeros([P cc]);      % means
Sw=zeros(P,P);          % Within class variability 


%-------------calculate Parameter-----------------
for i=1:cc;
    j = find(ctrain==i);                                     % select datapoints in this class
    n = length(j);                                           % number of sampels per category 
    muK(:,i) = sum(xtrain(:,j),2)/n;                         % get the Cluster means 
    res = bsxfun(@minus,xtrain(:,j),muK(:,i));
    Sw = Sw+res*res';                         % Estimate common covariance matrix
end;
%-------------Regularisation----------------------
P=size(Sw,1);
Sw=Sw/N; 
Sw=Sw+eye(P)*trace(Sw)*regularization/P;

%-------------classify----------------------------
[dummy N_xtest] = size(xtest);
L=muK'/Sw;                                  % Calculate the classifier L(i,:)=muK'*inv(Sw)
K=sum(L.*muK',2);                           % constant term for each class muK'*inv(Sw)*muK
G=bsxfun(@plus,-0.5*K,L*xtest);             % Classification function for each class, test point 
[gmax,idx]=max(G); 
cpred=idx;