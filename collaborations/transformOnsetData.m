function OnsetData = transformOnsetData(sn) %specify 1 through 24 for sn

baseDir= 'Z:/rhys/prepProd2/data'; %PC
saveDir= 'Z:/rhys/prepProd2/matlab/collaborations';

anaSubj = [3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42]; %Ps that reached performance threshold

subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24','s25','s26',...
    's27','s28','s29','s30','s31','s32','s33','s34','s35','s36','s37','s38','s39','s40','s41','s42'};
anaSn = anaSubj(sn);

cd(baseDir)
load(fullfile('imaging', 'GLM_firstlevel', subj_name{anaSn}, 'SPM.mat'), 'SPM')

%output must be trialN x 5 table containing:
%run       - run number
%cond      - regressor name (string)
%para      - press number for prod trials
%onsetSec  - onset in seconds relative to run start
%dur       - duration of regressor

OnsetData = []; %output variable, table

A = [];
behBlockNames = {'exp_BN47', 'exp_BN48', 'exp_BN49', 'exp_BN50', 'exp_BN51', 'exp_BN52'};

for i=1:length(behBlockNames)
    load(fullfile('behavioural', subj_name{anaSn}, [behBlockNames{i} '.mat']), 'B')
    A = addstruct(A, B);
    clear B
end

%%%For calculating trial length in error durations
A.goCueDur=nan(size(A.trialType,1),1);
A.goCueDur(A.trialType==1,1)=4000; %CHECK trial type!
A.goCueDur(A.trialType==2,1)=1000;

seqConditions = {...
    'O1T1 Prod', 'O1T2 Prod', 'O2T1 Prod', 'O2T2 Prod'; ...
    'O1T1 Prep', 'O1T2 Prep', 'O2T1 Prep', 'O2T2 Prep'...
    };

%%%Production trial regressors
prodOnsetSec = nan(size(A.timing));
for i=1:size(A.timing, 2)%for nPresses
    prodOnsetSec(:,i) = (A.tZero + A.timing(:,i)) / 1000;
end%for nPresses
prodRun  = kron(1:length(behBlockNames), ones(1,240))';
prodCond = cell(length(prodOnsetSec), 1);
prodCondIdx = nan(length(prodOnsetSec), 1);
prodCondIdx(~isnan(A.timing(:,1))) = A.seqID(~isnan(A.timing(:,1)));
for i=find(~isnan(A.timing(:,1)))'
    prodCond{i,1} = seqConditions(1, prodCondIdx(i));
end
prodCond     = repmat(prodCond, 1, 5);
prodCond     = reshape(prodCond', [], 1);
prodPara     = repmat(1:5, length(prodOnsetSec), 1);
prodPara     = reshape(prodPara', [], 1);
prodOnsetSec = reshape(prodOnsetSec', [], 1);
prodDur      = zeros(length(prodOnsetSec), 1);

%%%Preparation trial regressors
prepRun      = prodRun;
prepOnsetSec = nan(size(A.timing));
prepOnsetSec(isnan(A.timing(:,1))) = (A.tZero(isnan(A.timing(:,1)))) / 1000;
prepCond    = cell(length(prepOnsetSec), 1);
prepCondIdx = nan(length(prepOnsetSec), 1);
prepCondIdx(isnan(A.timing(:,1))) = A.seqID(isnan(A.timing(:,1)));
for i=find(isnan(A.timing(:,1)))'
    prepCond{i,1} = seqConditions(2, prepCondIdx(i));
end
prepCond     = repmat(prepCond, 1, 5);
prepCond     = reshape(prepCond', [], 1);
prepOnsetSec = reshape(prepOnsetSec', [], 1);
prepDur   = max(A.tZero- A.tZeroCue)/1000;
prepDur   = repmat(prepDur, size(prepOnsetSec));
prepPara = nan(length(prodOnsetSec), 1);

%%%Erroneous trial regressors
errorOnsetSec                = nan(size(A.timing));
errorOnsetSec(A.points == 0, 1) = (A.tZeroCue(A.points == 0))/1000;
errorDur                     = nan(size(A.timing));
errorDur(A.points == 0, 1)      = ( (A.tZero(A.points == 0)- A.tZeroCue(A.points == 0)) + A.goCueDur(A.points == 0) + A.crossDur(A.points == 0) + A.feedbackCalcDur(A.points == 0) + 1000 ) / 1000; %whole trial
errorOnsetSec = reshape(errorOnsetSec', [], 1);
errorDur      = reshape(errorDur', [], 1);
errorPara = nan(length(errorOnsetSec), 1);
errorCond     = repmat({'error'}, length(errorOnsetSec), 1);
errorRun      = prodRun;

%%%Preparation period in production trial regressors
prepInProdOnsetSec = nan(size(A.timing));
prepInProdOnsetSec(~isnan(A.timing(:,1))) = ( A.tZeroCue(~isnan(A.timing(:,1))) ) / 1000;
prepInProdOnsetSec = reshape(prepInProdOnsetSec', [], 1);
prepInProdDur  = nan(size(A.timing));
prepInProdDur(~isnan(A.timing(:,1)))   = ( A.tZero(~isnan(A.timing(:,1))) - A.tZeroCue(~isnan(A.timing(:,1))) ) / 1000;
prepInProdDur  = reshape(prepInProdDur', [], 1);
prepInProdRun  = prodRun;
prepInProdCond = repmat({'Prep in prod'}, length(prepInProdOnsetSec), 1);
prepInProdPara = prepPara;

%%%Feedback in every trial
feedbackOnsetSec      = nan(size(A.timing));
feedbackOnsetSec(:,1) = (A.tZero + A.goCueDur + A.crossDur + A.feedbackCalcDur) / 1000;
feedbackOnsetSec      = reshape(feedbackOnsetSec', [], 1);
feedbackDur      = nan(size(A.timing));
feedbackDur(:,1) = A.feedback;
feedbackDur      = reshape(feedbackDur', [], 1);
feedbackRun = prodRun;
feedbackCond = repmat({'feedback'}, length(feedbackOnsetSec), 1);
feedbackPara = nan(length(feedbackOnsetSec), 1);

%%%Concatenate all variables across regressors
run      = vertcat(prepRun, prodRun, prepInProdRun, errorRun, feedbackRun);                          % - run number
cond     = vertcat(prepCond, prodCond, prepInProdCond, errorCond, feedbackCond);                     % - regressor name (string)
para     = vertcat(prepPara, prodPara, prepInProdPara, errorPara, feedbackPara);                     % - press number for prod trials
onsetSec = vertcat(prepOnsetSec, prodOnsetSec, prepInProdOnsetSec, errorOnsetSec, feedbackOnsetSec); % - onset in seconds relative to run start
dur      = vertcat(prepDur, prodDur, prepInProdDur, errorDur, feedbackDur);                          % - duration of regressor

OnsetDataTemp = table(run, cond, para, onsetSec, dur);

%%%Sort by run and remove rows where onsetSec is NaN (from condition assignment earlier)
OnsetData = sortrows(OnsetDataTemp, 1); %sort by imaging run
OnsetData = rmmissing(OnsetData, 1, 'DataVariables', 'onsetSec');

for i=unique(OnsetData.run)' %sort onset times within runs
    OnsetData(OnsetData.run == i, 1:5) = sortrows(OnsetData(OnsetData.run == i, 1:5), 4);
end%for run

save(fullfile(saveDir, ['onsetData_' subj_name{anaSn}]), 'OnsetData')
cd(saveDir);
end
