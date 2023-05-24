function data=tempord3_combinedclass_corrected4Main(X,c,run,train,test);

% Remove the mean of the run for each
for i=1:max(run)
    j=find(run==i);
    % X(:,j)=X(:,j)-repmat(mean(X(:,j),2),1,length(j)); % SLOW
    X(:,j)=bsxfun(@minus,X(:,j),mean(X(:,j),2)); % MUCH FASTER
end;

[U,S,V]=svd(X,0); % Singular value decomposition
[data(1,1),cpred]=crossval_takemultipleout(@classify_lda_KclassesQuicker_tempord3_corrected4Main,S*V',c,train,test);
