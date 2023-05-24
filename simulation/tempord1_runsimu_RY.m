function [ACC]=tempord1_runsimu_RY(vtemp,vord,vinter,vnoise,n)
%%Simulates data and crossvalidates n times


% fprintf(['Runs simulation with: VTEMP=',num2str(vtemp),'  VORD=',num2str(vord),...
%     '  VINTER=',num2str(vinter),'  VNOISE=',num2str(vnoise),...
%     '  ITERATIONS=',num2str(n),'...\n']);  
D=[];

for i=1:n
    [X,Y,YtempMean,YcombiCorrected4Temp]=tempord1_makedata_RY('vtemp',vtemp,'vord',vord,'vinter',vinter,'vnoise',vnoise);
    [zAcc_ovrall, zAcc_ord, zAcc_temp, zAcc_int]=tempord1_classify_RY(X,Y,YtempMean,YcombiCorrected4Temp);
    
    %********* for dimensionality analysis:
    D.ACC(i,:)=zAcc_ovrall;
    %*********
    
    D.ACC(i,1)     = zAcc_ovrall; D.var(i,1)=1;
    D.ACC(i+n,1)   = zAcc_ord;    D.var(i+n,1)=2;
    D.ACC(i+n*2,1) = zAcc_temp;   D.var(i+n*2,1)=3;    
    D.ACC(i+n*3,1) = zAcc_int;    D.var(i+n*3,1)=4;
end

%%Display figure with barplot
% figure;
barplot(D.var,D.ACC,'split',D.var,'style_bold','leg',{'ovr', 'ord', 'temp' ,'inter'},'leglocation','north');
title(['# of iterations=',num2str(n),'  ORD=',num2str(vord),'  TEMP=',num2str(vtemp),'  INTER=',num2str(vinter),'  NOISE=',num2str(vnoise)]);
%axis([0 7 0 1.2]);   
ylabel('acc');
xlabel('factor');

ACC=D.ACC;

