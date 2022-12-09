function prepProd2_recodeS02

baseDir = 'D:\projects\rhys\prepProd2\data\behavioural\';
participant = 's02';

for i=0:52 %load blocks one by one
    fname=sprintf('exp_BN%02d.mat',(i));
    fileIn=fullfile(baseDir,participant,fname);
    load(fileIn);
    
    for j=1:length(B.TN) %for each block, re-calculate trial timings, presses, and points based on smoothed forces
        
        
        B.thresholdedSmoothedForces{j} = B.thresholdedForces{j}; %generate smoothed forces to save
        %         B.thresholdedSmoothedForces{j} = smoothdata(B.thresholdedForces(j),'gaussian',100); %generate smoothed forces to save
        
        %         threshold = newton2volt(2.5); %volt threshold for press calculated from Newtons
        threshold = 1.5292;
        
        %threshold guided by previous research (elife, Berlot)
        idx = B.thresholdedSmoothedForces{j} >= threshold; %index of 0s for baseline & 1s where our data goes above press threshold
        
        %% Calculate number of Presses
        numPress = 0; %create numPress variable to use in loop
        isPress = diff(idx); %find where our index increases from 0 to 1 per channel (press initiation)
        nChannels = 10;
        
        %%%% Loop through channels, counting number of times threshold is reached (number of presses)
        for k = 1:nChannels
            numPress = numPress + (sum(isPress(:,k)==1)); %record press initiation
        end
        n = numPress; %use as our 'n' variable for later point calculation
        
        %% Calculate Press Timings
        
        %create five variables containing single channel press initiation time data, in channel order...
        %using our 'isPress' variable from before as a reference to choose timestamps
        thumbTemp = B.timeStamps{j}(isPress(:,1)==1)';
        thumb = [thumbTemp*1000; repelem(1,numel(thumbTemp))];
        
        indexTemp = B.timeStamps{j}(isPress(:,2)==1)';
        index = [indexTemp*1000; repelem(2,numel(indexTemp))];
        
        middleTemp = B.timeStamps{j}(isPress(:,3)==1)';
        middle = [middleTemp*1000; repelem(3,numel(middleTemp))];
        
        ringTemp = B.timeStamps{j}(isPress(:,4)==1)';
        ring = [ringTemp*1000; repelem(4,numel(ringTemp))];
        
        pinkyTemp = B.timeStamps{j}(isPress(:,5)==1)';
        pinky = [pinkyTemp*1000; repelem(5,numel(pinkyTemp))];
        
        presses = ([thumb index middle ring pinky]);% concatenate all first press times...
        presses(1,:) = presses(1,:)  - (B.transducerOffsetMeasured{j}); %and subtract delay caused by transducer starting earlier than go cue
        presses = sortrows(presses.',1).'; %sort our presses by timing (column-wise)
        
        keytime = presses(1,:);
        key = presses(2,:);
        
        %if no presses, fill the vector with NaNs instead
        if isempty(key)
            key=nan(1,5);
        end
        
        finger=key;
        
        for k=1:length(key)
            if isempty(find(key(k)==exp.key.map, 1)) % No response
                finger(k)=NaN;
                keytime(k)=NaN;
            else
                finger(k)=find(key(k)==exp.key.map);
            end
        end
        
        key=finger;
        
        
        
        % Probe trial reaction time
        catchRT=keytime(1);
        % Check actual delay period
        cueDurActual=tZero-tZeroCue; % Updated: duration between Go cue and fractal cue display
        
        %% Log responses (press, npress, timing)
        
        exp.sequence.press(trials,:)=nan(1,size(timing,2));
        exp.sequence.timing(trials,:)=nan(1,size(timing,2));
        
        exp.sequence.npress(trials,:)=n;
        if n==size(timing,2) %if correct number of responses
            exp.sequence.press(trials,:)=key;
            exp.sequence.timing(trials,:)=keytime';
            responseSwitch=1;
        elseif n>size(timing,2)
            exp.sequence.press(trials,:)=key(1:size(timing,2)); % Log the first five responses
            exp.sequence.timing(trials,:)=keytime(1:size(timing,2));
            responseSwitch=0;
        elseif n<size(timing,2)
            exp.sequence.press(trials,1:length(key))=key; % Log the first five responses
            exp.sequence.timing(trials,1:length(key))=keytime';
            responseSwitch=0;
        end
        
        %% Calculate Performance
        
        seqPerformance=sequence-exp.sequence.press(trials,:); %
        correctSum=nansum(seqPerformance==0);
        seqPerformance=2-(seqPerformance==0); %2 - wrong; 1 - correct;
        feedbackSym={'x','-'};
        for j=1:length(seqPerformance)
            seqPerformanceFeedback(j)=feedbackSym{seqPerformance(j)};
        end
        seqPerformanceFeedbackString=strjoin(seqPerformanceFeedback,' ');
        seqPerformance=(1-nansum(seqPerformance~=0)/length(sequence)); %every error substracted from 4
        
        %% Points and Feedback %%%%%%%%%%%%%%%%%%%%%%
        if T.trialType(trials,1)==1 %%% Sequence trials
            disp('Sequence');
            timingPerformance=exp.sequence.timing(trials,:)-exp.sequence.cueTime(trials,:);
            timingPerformanceRT=exp.sequence.timing(trials,1); %seq RT for 1st press
            timingPerformance=timingPerformance(2:5) ./ diff([0 exp.sequence.cueTime(trials,2:5)]) * 100; % Deviation from target timing as percent of target interval
            timingPerformance=[0 timingPerformance];
            %disp(['timingPerformance: ' num2str(timingPerformance)]);
            %%%% Visual FB
            limLines=['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'];% display limits
            centreLines=['-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'];
            preparestring(limLines,(i+3),0, -100);
            preparestring(limLines,(i+3),0, +100);
            preparestring(centreLines,(i+3),0, 0);
            
            %%%% Points for sequence RT:
            targetSum=5;
            if timingPerformanceRT>=0 && timingPerformanceRT<=200 && correctSum==targetSum
                exp.sequence.points(trials,1)=5;
            elseif timingPerformanceRT>200 && timingPerformanceRT<=360 && correctSum==targetSum
                exp.sequence.points(trials,1)=4;
            elseif timingPerformanceRT>360 && timingPerformanceRT<=480 && correctSum==targetSum
                exp.sequence.points(trials,1)=3;
            elseif timingPerformanceRT>480 && timingPerformanceRT<=560 && correctSum==targetSum
                exp.sequence.points(trials,1)=2;
            elseif timingPerformanceRT>560 && timingPerformanceRT<=600 && correctSum==targetSum
                exp.sequence.points(trials,1)=1;
            elseif timingPerformanceRT>600 && correctSum==targetSum
                exp.sequence.points(trials,1)=0;
            elseif timingPerformanceRT<0 && correctSum==targetSum % if they press a button prematurely before the digit cue (tZero)
                exp.sequence.points(trials,1)=0;
            elseif correctSum~=targetSum
                exp.sequence.points(trials,1)=0;
            end
            seqRTpoints=exp.sequence.points(trials,1);
            disp(['RT Points: ' num2str(seqRTpoints)]);
            
            %%%% Points for percent deviation from target interval:
            timingPerformanceInt=timingPerformance(2:5);
            if responseSwitch==0 % Wrong number of responses
                exp.sequence.points(trials,1)=0;
                seqRTpoints=0;
            elseif nanmean(abs(timingPerformanceInt))>=0 && nanmean(abs(timingPerformanceInt))<=10 && correctSum==targetSum
                exp.sequence.points(trials,1)=5;
            elseif nanmean(abs(timingPerformanceInt))>10 && nanmean(abs(timingPerformanceInt))<=20 && correctSum==targetSum
                exp.sequence.points(trials,1)=4;
            elseif nanmean(abs(timingPerformanceInt))>20 && nanmean(abs(timingPerformanceInt))<=30 && correctSum==targetSum
                exp.sequence.points(trials,1)=3;
            elseif nanmean(abs(timingPerformanceInt))>30 && nanmean(abs(timingPerformanceInt))<=40 && correctSum==targetSum
                exp.sequence.points(trials,1)=2;
            elseif nanmean(abs(timingPerformanceInt))>40 && nanmean(abs(timingPerformanceInt))<=50 && correctSum==targetSum
                exp.sequence.points(trials,1)=1;
            elseif nanmean(abs(timingPerformanceInt))>50 && correctSum==targetSum
                exp.sequence.points(trials,1)=0;
            elseif correctSum~=targetSum % Incorrect responses
                exp.sequence.points(trials,1)=0;
            end
            
            if exp.sequence.timing(trials,1)<0 % if they press a button prematurely before the Go cue
                exp.sequence.points(trials,1)=0;
            end
            
            seqIntpoints=exp.sequence.points(trials,1);
            disp(['Int Points: ' num2str(seqIntpoints)]);
            seqTotalpoints=seqIntpoints+seqRTpoints; % Per trial
            % Map points to auditory feedback
            if seqTotalpoints==0
                soundName=soundFeedback(1);
            elseif seqTotalpoints==1
                soundName=soundFeedback(2);
            elseif seqTotalpoints==2
                soundName=soundFeedback(3);
            elseif seqTotalpoints==3
                soundName=soundFeedback(4);
            elseif seqTotalpoints==4
                soundName=soundFeedback(5);
            elseif seqTotalpoints==5
                soundName=soundFeedback(6);
            elseif seqTotalpoints==6
                soundName=soundFeedback(7);
            elseif seqTotalpoints==7
                soundName=soundFeedback(8);
            elseif seqTotalpoints==8
                soundName=soundFeedback(9);
            elseif seqTotalpoints==9
                soundName=soundFeedback(10);
            elseif seqTotalpoints==10
                soundName=soundFeedback(11);
            end
            disp(['Points: ' num2str(seqTotalpoints)]);
            pointsAll=pointsAll+seqTotalpoints; % Update points
            exp.sequence.points(trials,1)=seqIntpoints+seqRTpoints; % Log points
            
        else   %%% catch trials
            disp('Catch');
            catchRT=keytime(1);
            targetSum=0;
            % Reaction time (ms):
            if n>0 % any response recorded
                exp.sequence.points(trials,1)=0;
                soundName=soundFeedback(1);
            elseif n==targetSum % no response
                exp.sequence.points(trials,1)=5;
                soundName=soundFeedback(6);
            elseif catchRT<0 % if they press a button prematurely before the noGo cue (tZero)
                exp.sequence.points(trials,1)=0;
                soundName=soundFeedback(1);
            elseif correctSum~=targetSum
                exp.sequence.points(trials,1)=0;
                soundName=soundFeedback(1);
            end
            catchRTpoints=exp.sequence.points(trials,1);
            disp(['Points: ' num2str(catchRTpoints)]);
            pointsAll=pointsAll+catchRTpoints; % Update points
            exp.sequence.points(trials,1)=catchRTpoints; % Log points
        end
        
    end
    
    
end

disp('all done')