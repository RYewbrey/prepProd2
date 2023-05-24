function data=tempord3_combinedclass(X,c,run,train,test,varargin);

% Remove the mean of the run for each
for i=1:max(run)
    j=find(run==i);
    % X(:,j)=X(:,j)-repmat(mean(X(:,j),2),1,length(j)); % SLOW
    X(:,j)=bsxfun(@minus,X(:,j),mean(X(:,j),2)); % MUCH FASTER
end;


%%Mean searchlight! MVPA only on the mean across searchlight - no patterns;

%%

[U,S,V]=svd(X,0); % Singular value decomposition

[data(1,1),cpred]=crossval_takemultipleout(@classify_lda_KclassesQuicker,S*V',c,train,test);
% data=repmat(data,3,1);
        
%%% Dimensionality analysis
%  dim=varargin{1};
% [data(1:8),cpred]=crossval_takemultipleout_multi(@classify_lda_Kclasses_dim,S*V',c,train,test,[1:8]);
%  data=data(dim);




%Substract out the mean;
%Comment out for whole brain analysis as in Kornysheva&Diedrichsen 2013
%                  [data(1,1),cpred]=crossval_takemultipleout(@classify_lda_KclassesQuicker_withoutMean,X,c,train,test);
                
  
%Mean searchlight! MVPA only on the mean across searchlight - no patterns;
%Comment out for whole brain analysis as in Kornysheva&Diedrichsen 2013
%              [data(1,1),cpred]=crossval_takemultipleout(@classify_lda_KclassesQuicker_Mean,S*V',c,train,test); 
%                 data=repmat(data,3,1);


% cTemp=floor((c-1)/3)+1;
% cOrd=mod((c-1),3)+1;
% cTempP=floor((cpred-1)/3)+1;
% cOrdP=mod((cpred-1),3)+1;
% data(2,1)=sum(cTemp==cTempP)/length(cTemp); % did it get the temporal structure right
% data(3,1)=sum(cOrd==cOrdP)/length(cOrd);   % did it get the ordinal structure right