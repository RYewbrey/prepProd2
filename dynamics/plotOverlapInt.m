function plotOverlapInt(trial)

fprintf('Trial %d\n', trial)
evalin('base', 'save tempvars.mat') %works on current workspace variables *to change*

load tempvars.mat


peakChannels = [1 2; 3 4; 5 6; 7 8]; %channels contained within each peak as a row

figure
subplot(2, 4, [1 4]), plot(F.threshForce{trial}) %plot force traces
plotYLim = axis;
plotYLim = plotYLim(4); %variable to hold the range of the y axis of force trace plots

for i = 1:4 %loop through all intervals of a trial
    
    fprintf('Overlap on interval %d = %d\n',i, F.areaInt(trial, i))
    
    subplot(2, 4, i+4), plot(F.peaks{trial, peakChannels(i, 1)}, 'b')
    ylim([0 plotYLim]) %plot and set limit of y axis to equal of force traces
    
    hold on
    subplot(2, 4, i+4), plot(F.peaks{trial, peakChannels(i, 2)}, 'g')
    ylim([0 plotYLim])
    subplot(2, 4, i+4), plot(F.overlapMeasure{trial, i}, 'r')
    ylim([0 plotYLim])
end