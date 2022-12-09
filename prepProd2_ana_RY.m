function D=prepProd2_ana_RY(sub,baseDir,varargin)
% Experiment prepProd - SAMlab©2019

%%Script to produce all of the variables required for prepProd_ana_run%%

% Analysis of subject RTs and accuracies on the spatio-temporal production task in the behavioural training

subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

%%% Blocks:
% BN01-BN04 probePre blocks
% BN05-BN22 training Day1
% BN23-BN40 training Day2
% BN29-BN40 Production from memory only (performance check!)
% BN41-BN44 probePost blocks
% BN45-BN46 training refresher
% BN47-BN52 fMRI

%%% trialType coding:
% trialType==1 sequence (instructed & mem)
% trialType==2 catch
%%% mode coding:
% mode==1 sequence instructed
% mode==2 sequence mem
% mode==0 catch

%% Blocks to Analyse
minBN=1;
maxBN=52;
% maxBN=52;
% vararginoptions(varargin,{'minBN','maxBN'});

% analysis of the spatial and temporal accuracy of subject performance 

D=[];
for i=minBN:maxBN
    fname=sprintf('exp_BN%02d.mat',(i));
    fileIn=fullfile(baseDir,subj_name{sub},fname);
    load(fileIn);
    
    B.timingJitter=abs(max(B.cueTime-B.cueTimePlanned,[],2));%std(B.cueTime-B.cueTimePlanned,[],2);
    B.timingJitterInterval=mean(abs(diff([zeros(size(B.cueTime,1),1) B.cueTime],1,2)-diff([zeros(size(B.cueTimePlanned,1),1) B.cueTimePlanned],1,2)),2);
    B.jitterMask=zeros(size(B.timingJitter,1),1);
    B.jitterMask(B.timingJitterInterval>10)=1; %more than 10ms difference to target interval
%     B.timingInterval=diff([zeros(size(B.cueTimePlanned,1),1) B.timing],1,2);   %old
    B.timingIntervalTarget=diff([zeros(size(B.cueTimePlanned,1),1) B.cueTimePlanned],1,2);
    B.RT=nanmean((B.timing-B.cueTimePlanned),2);
    B.RT(B.trialType==2)=B.timing(B.trialType==2,1); %catch trials in case of a recorded response
    B.RTabs=nanmean(abs(B.timing-B.cueTimePlanned),2); %pre vs post test - mean RT deviation in ms 
    B.RTabsFinger=abs(B.timing-B.cueTimePlanned); % RT deviation from target in ms for each press 
    % Percent deviation RT from target interval
    B.RTabsInterval=nanmean(abs((diff([B.timing(:,1) B.timing(:,2:5)],1,2)-diff([B.cueTimePlanned(:,1) B.cueTimePlanned(:,2:5)],1,2))) ./ diff([B.cueTimePlanned(:,1) B.cueTimePlanned(:,2:5)],1,2),2)*100; %corrected
    B.RTabsFingerInterval=abs((diff([B.timing(:,1) B.timing(:,2:5)],1,2)-diff([B.cueTimePlanned(:,1) B.cueTimePlanned(:,2:5)],1,2))) ./ diff([B.cueTimePlanned(:,1) B.cueTimePlanned(:,2:5)],1,2)*100; %updated Katja

    %% intervals between finger presses
    B.timingInterval=diff(B.timing,1,2); %interval to Go cue is not interesting
    B.timingInterval1=B.timingInterval(:,1);
    B.timingInterval2=B.timingInterval(:,2);
    B.timingInterval3=B.timingInterval(:,3);
    B.timingInterval4=B.timingInterval(:,4);
    
    B.RTabsFingerBegin=(abs(B.timing-B.cueTimePlanned) ./ B.cueTimePlanned) * 100;
    %     B.RTFingerBegin=((B.timing-B.cueTimePlanned) ./ B.cueTimePlanned) * 100+100;
    
    B.targetIntervalDiff=diff([B.cueTimePlanned(:,1) B.cueTimePlanned(:,2:5)],1,2);
    B.targetIntervalDiff1=B.targetIntervalDiff(:,1);
    B.targetIntervalDiff2=B.targetIntervalDiff(:,2);
    B.targetIntervalDiff3=B.targetIntervalDiff(:,3);
    B.targetIntervalDiff4=B.targetIntervalDiff(:,4);
    
    D=addstruct(D,B,'row');
    B=[];
end;

D.pointsFinal=nan(size(D.timing,1),1);
D.pointsFinal=max(D.pointsAll,1); % take max for each subj

%% Error:incorrect presses
errorFinger=(D.cueFinger-D.press);
errorFinger(isnan(errorFinger))=1; % true (1) where NaN
D.errorFinger=any(errorFinger,2); % true (1) where incorrect press (nonzero)
D.errorFinger(D.trialType==2,:)=errorFinger(D.trialType==2,1);

%% Accuracy for probes
D.correctFinger=~D.errorFinger; % invert logical true false vector
D.correctFinger(D.trialType==2,:)=D.correctFinger(D.trialType==2,1);

%% Accuracy for memory sequences
errorFingerSeq=(D.cueFinger-D.press);
errorFingerSeq(isnan(errorFingerSeq))=1; % true (1) where NaN
for i=1:5 % for each column
y=errorFingerSeq(:,i); 
D.errorFingerSeq(:,i)=any(y,2); % assign true (1) where incorrect press on each col
end
D.errorFingerSeq(D.trialType==1&D.mode==2,:)=errorFingerSeq(D.trialType==1&D.mode==2,:); % 
D.correctFingerSeq=~D.errorFingerSeq; % invert logical true false vector for accuracy (true 1 where correct)
D.correctFingerSeq(D.trialType==1&D.mode==2,:)=D.correctFingerSeq(D.trialType==1&D.mode==2,:);
% D.errorFingerSeq=errorFingerSeq(:); % convert matrix to a vector

%% Error for memory sequences
errorFingerSeq=(D.cueFinger-D.press);
errorFingerSeq(isnan(errorFingerSeq))=1; % true (1) where NaN
for i=1:5 % for each column
y=errorFingerSeq(:,i); 
D.errorFingerSeq(:,i)=any(y,2); % assign true (1) where incorrect press on each col
end
D.errorFingerSeq(D.trialType==1&D.mode==2,:)=errorFingerSeq(D.trialType==1&D.mode==2,:); % 

%% Conditions:
%Trained vs. controls
   % sequence conditions pooled here
D.seqCond=nan(size(D.BN,1),1);
D.seqCond(D.seqID>=1&D.seqID<=4,1)=1; % trained seq (green)
D.seqCond(D.seqID>=5&D.seqID<=8,1)=2; % control temp transfer (red)
D.seqCond(D.seqID>=9&D.seqID<=12,1)=3; % control spatial transfer (blue)
D.seqCond(D.seqID>=13&D.seqID<=16,1)=4; % control new(black)

%Pre vs. post
D.prepost=nan(size(D.BN,1),1);
D.prepost(D.BN>=1&D.BN<=4,:)=1; % pre training
D.prepost(D.BN>=41&D.BN<=44,:)=2; % post

%Training Day1 vs. Day2
D.training=nan(size(D.BN,1),1);
D.training(D.BN>=5&D.BN<=22,:)=1; % Day1
D.training(D.BN>=23&D.BN<=40,:)=2; % Day2

%Sequence cue + fix cross preparation duration
D.durCond=nan(size(D.cueDur,1),1);
D.durCond(D.cueDur==983,1)=1; 
D.durCond(D.cueDur==1500,1)=2;
D.durCond(D.cueDur==1983,1)=3;
D.durCond(D.cueDur==2500,1)=4;

%%% Mean delay time for each delay condition
cueDurT1mean=round(mean(D.cueDurActual(D.durCond==1,:))); % get the mean for each
cueDurT2mean=round(mean(D.cueDurActual(D.durCond==2,:)));
cueDurT3mean=round(mean(D.cueDurActual(D.durCond==3,:)));
cueDurT4mean=round(mean(D.cueDurActual(D.durCond==4,:)));
D.cueDurActual(D.durCond==1,:)=cueDurT1mean;
D.cueDurActual(D.durCond==2,:)=cueDurT2mean;
D.cueDurActual(D.durCond==3,:)=cueDurT3mean;
D.cueDurActual(D.durCond==4,:)=cueDurT4mean;
D.durCond(D.cueDurActual==cueDurT1mean,1)=1; % use the actual mean duration of delay for each condition
D.durCond(D.cueDurActual==cueDurT2mean,1)=2;
D.durCond(D.cueDurActual==cueDurT3mean,1)=3;
D.durCond(D.cueDurActual==cueDurT4mean,1)=4;

        %% Each IPI expressed as percent of the mean long duration
%%% Mean of reference condition across intervals for each subj
intervalRef=D.timingInterval(D.trialType==1&D.mode==2&D.durCond==3&D.BN>=47,:); %take four IPIs 
meanIntervalRef=nanmean(intervalRef(:)); %collate columns into a single vector and take the mean
% targetIntervalRef=repmat([100 100 100],length(intervalRef),1);
targetIntervalRef=D.targetIntervalDiff(D.trialType==1&D.mode==2&D.durCond==3&D.BN>=47,:);
meanTargetIntervalRef=nanmean(targetIntervalRef(:));
D.intervalPrc=nan(size(D.timingInterval,1),4);
D.intervalPrc1=nan(size(D.timingInterval,1),1);
D.intervalPrc2=nan(size(D.timingInterval,1),1);
D.intervalPrc3=nan(size(D.timingInterval,1),1);
D.intervalPrc4=nan(size(D.timingInterval,1),1);
D.targetIntervalPrc=nan(size(D.timingInterval,1),4);
D.targetIntervalPrc1=nan(size(D.timingInterval,1),1);
D.targetIntervalPrc2=nan(size(D.timingInterval,1),1);
D.targetIntervalPrc3=nan(size(D.timingInterval,1),1);
D.targetIntervalPrc4=nan(size(D.timingInterval,1),1);
D.timErrorRelAbs=nan(size(D.timingInterval,1),4);
D.timErrorRel=nan(size(D.timing,1),1);
%For each condition calculate each actual and target IPI relative to the above reference mean & get the percent:
for i=1:4 %delay
    for j=1:4 %interval/column
    D.intervalPrc(D.trialType==1&D.mode==2&D.durCond==i&D.BN>=47,j)=D.timingInterval(D.trialType==1&D.mode==2&D.durCond==i&D.BN>=47,j)/meanIntervalRef*100;
    D.targetIntervalPrc(D.trialType==1&D.mode==2&D.durCond==i&D.BN>=47,j)=D.targetIntervalDiff(D.trialType==1&D.mode==2&D.durCond==i&D.BN>=47,j)/meanTargetIntervalRef*100;
    end
end
%Next,calculate deviation of interval prc from target interval prc for each trial:
D.timErrorRelAbs(D.trialType==1&D.mode==2&D.BN>=47,:)=abs(D.intervalPrc(D.trialType==1&D.mode==2&D.BN>=47,:) - D.targetIntervalPrc(D.trialType==1&D.mode==2&D.BN>=47,:)); %get absolute difference
D.timErrorRel(D.trialType==1&D.mode==2&D.BN>=47,1)=nanmean(D.timErrorRelAbs(D.trialType==1&D.mode==2&D.BN>=47,:),2); %get the mean per trial/row
D.intervalPrc1(D.trialType==1&D.mode==2&D.BN>=47,1)=D.intervalPrc(D.trialType==1&D.mode==2&D.BN>=47,1);
D.intervalPrc2(D.trialType==1&D.mode==2&D.BN>=47,1)=D.intervalPrc(D.trialType==1&D.mode==2&D.BN>=47,2);
D.intervalPrc3(D.trialType==1&D.mode==2&D.BN>=47,1)=D.intervalPrc(D.trialType==1&D.mode==2&D.BN>=47,3);
D.intervalPrc4(D.trialType==1&D.mode==2&D.BN>=47,1)=D.intervalPrc(D.trialType==1&D.mode==2&D.BN>=47,4);
D.targetIntervalPrc1(D.trialType==1&D.mode==2&D.BN>=47,1)=D.targetIntervalPrc(D.trialType==1&D.mode==2&D.BN>=47,1);
D.targetIntervalPrc2(D.trialType==1&D.mode==2&D.BN>=47,1)=D.targetIntervalPrc(D.trialType==1&D.mode==2&D.BN>=47,2);
D.targetIntervalPrc3(D.trialType==1&D.mode==2&D.BN>=47,1)=D.targetIntervalPrc(D.trialType==1&D.mode==2&D.BN>=47,3);    
D.targetIntervalPrc4(D.trialType==1&D.mode==2&D.BN>=47,1)=D.targetIntervalPrc(D.trialType==1&D.mode==2&D.BN>=47,4);

%%% Create vector of 8 repetitions per seqID for plotting RTabs(ms) as a function of repetition
maxRep=8; % number of successive repetitions for each seqID 
D.trialRep=repmat([1:maxRep]',size(D.TN,1)/maxRep,1);

%Calculate error rate for prepost-seqCond
D.errorRate=nan(size(D.BN,1),1);
for j=1:max(D.prepost)
    for k=1:max(D.seqCond)
        D.errorRate(D.prepost==j&D.seqCond==k,1)=(sum(D.errorFinger(D.prepost==j&D.seqCond==k),1))/size(D.errorRate(D.prepost==j&D.seqCond==k),1)*100 ;
    end;
end;
%Calculate error rate for training phase
for j=1:max(D.training)
    for k=1:max(D.seqCond)
        D.errorRate(D.training==j&D.seqCond==k,1)=(sum(D.errorFinger(D.training==j&D.seqCond==k),1))/size(D.errorRate(D.training==j&D.seqCond==k),1)*100 ;
    end;
end

D.TNall=(1:size(D.RT,1))';
% 
% D.cumPoints=cumsum(D.points);
% 
% 
% D;
