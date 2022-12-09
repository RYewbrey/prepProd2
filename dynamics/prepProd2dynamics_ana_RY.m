function prepProd2dynamics_ana_RY(subj)

%Quantify amount of overlap between consecutive presses: Core script (generates key variables)
%
%PERFORMS:
%
%   Loop through blocks, trials.
%       Spit out area under the curve calculations within participants, specified by varargin.
%   Determine finger order from cueFinger in behavioural file.
%   Saves Force thresholded force data and peak calculations.
%
%--------------------------------------------------------------------------------------------------------------------------

%%% Add before starting any scripts (comment in when pasting into command line):
% addpath(genpath('E:\projects\rhys\prepProd2\matlab\dynamics'));
% addpath(genpath('E:\projects\rhys\prepProd2\data\behavioural\forces')); %path to behavioural data containing force traces
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)

baseDir = 'E:\projects\rhys\prepProd2\data\behavioural\forces'; %save location of force data
groupDir = 'E:\projects\rhys\prepProd2\data\behavioural\forces\group'; %where to save processed data

if subj == 1 || subj == 2 || subj == 8 || subj == 14 || subj == 19 || subj == 23 || subj == 24 || subj == 27 || subj == 28 || subj == 29 || subj == 30 || subj == 35; %subjects who did not reach final session
    disp('Invalid subject number. Please check subjects to be analysed.')
    return %leave script if an invalid subj is selected
end

subjName={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

blockN = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', ... %block numbers
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', ...
    '26', '27','28','29','30','31','32','33','34','35','36','37','38','39',  ...
    '40','41','42','43','44','45','46','47','48','49','50','51','52'};

cd([baseDir, '\', subjName{subj}])

for i=1:length(blockN) %loop through subject's blocks and load into workspace 'forceData' variable
    datafilename=['exp_BN' blockN(i) '.mat'];
    datafilename = strjoin(datafilename, ''); %converts from cell to string
    
    try
        if i == 1 %on first loop, create 'move' variable
            load(datafilename) %load block
            forceData = E;
        else
            load(datafilename)
            forceData = addstruct(forceData, E);
        end
    catch
        disp(sprintf('%s not loaded. Check file location and presence.', datafilename))
    end
end

if (isempty(forceData)) %if move contains nothing, exit script
    disp('No data available')
    return;
end;

F = forceData; %create F struct to contain subject-specific analysis

% By this point, F is a struct containing all variables for all of a participant's trials

for trial=1:length(F.thresholdedSmoothedForces)%for every trial, we want a value of overlap for each consecutive press 
    
    order = F.cueFinger(trial, :); %determine finger order from .mat file
        
    prodChannels = 1:5; %force channels corresponding to production fingers
    
    force = F.thresholdedSmoothedForces{trial}(:, prodChannels); %load channels of target trial into new force matrix
    
    if isnan(order) %for no-go trials, assign nan to all variables
        
        R.orderedThresholdedForce{trial,:} = NaN;
        
        data = cell(1,8); %generate 1x8 NaN cells to fill peaks
        data(1:8) = {NaN(1,1)};
        R.peaks(trial,:) = data;
        
        R.overlapMeasure (trial,:) = data(1:4);
        R.areaInt(trial,:) = NaN(1,4);
        R.areaIntAvg(trial) = NaN;
        
    else %for production trials, run overlap calcs
        
        %thresholding
        force1 = force(:,order(1)); %new variables for each channel, ordered according to ordseq
        force2 = force(:,order(2));
        force3 = force(:,order(3));
        force4 = force(:,order(4));
        force5 = force(:,order(5));
        
        orderedForce = [force1 force2 force3 force4 force5];
        orderedForce(orderedForce < 0) = 0; %set lowest value to 0, remove negative readings
        
        targetChannels = [1 2; 2 3; 3 4; 4 5];
        peakCounter = 1;
        
        for intervals = 1:4
            
            peak = zeros(1,2); %pre-allocate variable for loop
            j=1; %within-loop counter
            for i = targetChannels(intervals,:) %find the peaks of two adjacent presses
                
                [~, peak(j)] = max(orderedForce(:,i)); %find the location of the press peaks and assign them to 'peak' variable
                j=j+1;
                
            end
            
            R.orderedThresholdedForce{trial,:} = orderedForce;
            
            threshPeak1 = orderedForce(:, targetChannels(intervals,1));
            threshPeak2 = orderedForce(:, targetChannels(intervals,2));
            
            peak1 = threshPeak1(peak(1):peak(2))';
            peak2 = threshPeak2(peak(1):peak(2))';
            
            R.peaks{trial, peakCounter} = peak1;
            peakCounter = peakCounter + 1;
            
            R.peaks{trial, peakCounter} = peak2;
            peakCounter = peakCounter + 1;
            
            overlapPlot = [peak2(peak2<peak1) peak2(peak2==peak1) peak1(peak1<peak2)];
            R.overlapMeasure{trial, intervals} = overlapPlot;
            
            R.areaInt(trial,intervals) = trapz(overlapPlot);
        end
        
        R.areaIntAvg(trial,:) = mean(R.areaInt(trial,:));        
    end
end

%--------------------- Save variables of interest in their whole form, outside the loop ---------------------

%% Error:incorrect presses
errorFinger=(F.cueFinger-F.press);
errorFinger(isnan(errorFinger))=1; % true (1) where NaN
F.errorFinger=any(errorFinger,2); % true (1) where incorrect press (nonzero)
F.errorFinger(F.trialType==2,:)=errorFinger(F.trialType==2,1);

R.errorFinger = F.errorFinger; %1 = erroneous trial, 0 = correct trial

R.mode = F.mode; %0 = no-production, 1 = instructed, 2 = non-instructed
R.targetPressOrder = F.cueFinger;
R.actualPressOrder = F.press;
R.thresholdedSmoothedForces = F.thresholdedSmoothedForces;
R.seqID = F.seqID;
R.BN = F.BN;
R.subj = F.subj;

if ~exist(groupDir, 'dir')
       mkdir(groupDir)
end

cd(groupDir)
disp('Saving...')
save(sprintf('%s_overlap_data', subjName{subj}),'R') %area_int, peaks, overlapPlot, thresholded force data to .mat file
disp(sprintf('%s done', subjName{subj}))

%--------------------------------------------------------------------------------------------------------------------------

%% Example of function

% % Sample Data
% x1 = 0.01:0.01:1;
% y1 = [sin(pi*x1) zeros(1,20)];
% x2 = x1+0.2;
% y2 = [zeros(1,20) sin(pi*x1)];
% y_d = [y2(y2<y1) y1(y1<y2)];
% area_int  = trapz(y_d)
% %plots in case you want to visualize
% plot(y1, 'red')
% hold on
% plot(y2, 'blue')
% plot(y_d,'k-o')



% %Plot each finger from trial seperately
% figure
%
% % numFingers = 10;
% counter = 0;
%
% for i = 4:13
%     counter = counter + 1;
%     subplot(5, 5, counter), plot(MOV{1,1} (:,i))
%
%     axis([0 1302 0 2])
% end