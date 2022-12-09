function prepProdSimu_runsimuPP(prepProd,vtemp,vord,vinter,vnoise,n)
%==============================================================================
% Simulates data across n iterations
%
% [ACC] = prepProdSimu_runsimu(vtemp,vord,vinter,vnoise,n)
% prepProd = 0: prep=prod, 1: prep~=prod
% vtemp = variance across timing conditions - higher variance = better decoding
% vord = variance across order conditions
% vinter = variance across four sequences, for overall decoding
% n = number of iterations
%
% Typical variance from 0.01 - 1, but can increase above that. Check Kornysheva
% & Diedrichsen elife 2014 for reasonable decoder accuracy values
% RY 09/2022

% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('D:\projects\toolboxes\pcm_toolbox')); %PCM toolbox

saveDir = 'G:\projectsBackup\rhys\prepProd2\data\imaging\simulations';

D=[];
phase = {'prep', 'prod'};
nseq = 4; %number of sequences collapsed across phase

% for p=1:length(phase)%for phases
for i=1:n %for iterations
    Y=prepProdSimu_makedataPP('prepProd',prepProd,'vtemp',vtemp,'vord',vord,'vinter',vinter,'vnoise',vnoise);
    
    prepData = [Y(1:4,:); Y(9:12,:); Y(17:20,:); Y(25:28,:); Y(33:36,:); Y(41:44,:)];
    prodData = [Y(5:8,:); Y(13:16,:); Y(21:24,:); Y(29:32,:); Y(37:40,:); Y(45:48,:)];
    
    [acc_ovrall_prep,acc_int_prep, acc_ord_prep, acc_temp_prep]=prepProdSimu_classify(prepData);
    [acc_ovrall_prod,acc_int_prod, acc_ord_prod, acc_temp_prod]=prepProdSimu_classify(prodData);
    
    D.acc(i,1) = acc_ovrall_prep; D.var(i,1)=1; D.phase(i,1)=1;
    D.acc(i,1) = acc_ovrall_prod; D.var(i,1)=1; D.phase(i,1)=2;
    
    D.acc(i+n,1) = acc_int_prep; D.var(i+n,1)=2; D.phase(i+n,1)=1;
    D.acc(i+n,1) = acc_int_prod; D.var(i+n,1)=2; D.phase(i+n,1)=2;
    
    D.acc(i+n*2,1)= acc_ord_prep; D.var(i+n*2,1)=3; D.phase(i+n*2,1)=1;
    D.acc(i+n*2,1)= acc_ord_prod; D.var(i+n*2,1)=3; D.phase(i+n*2,1)=2;
    
    D.acc(i+n*3,1)= acc_temp_prep; D.var(i+n*3,1)=4; D.phase(i+n*3,1)=1;
    D.acc(i+n*3,1)= acc_temp_prod; D.var(i+n*3,1)=4; D.phase(i+n*3,1)=2;
    
    D.data{i,1} = Y;
    
end%for iterations
% end%for phases
%Overall and integrated z value transformation
numTests=6;
numCat=4;
mu=1/numCat; %mu=0.25;
N=numTests*numCat;
sigma=sqrt(mu*(1-mu)*1/N);

D.acc(D.var == 1) = (D.acc(D.var == 1)-mu)/sigma;
D.acc(D.var == 2) = (D.acc(D.var == 2)-mu)/sigma;

%order and timing z value transformation
takeOneOutIter=2;
numTests=1;
numCat=2;
mu=1/numCat; %mu=0.5;
N=numTests*numCat*takeOneOutIter;
sigma=sqrt(mu*(1-mu)*1/N);

D.acc(D.var == 3) = (D.acc(D.var == 3)-mu)/sigma;
D.acc(D.var == 4) = (D.acc(D.var == 4)-mu)/sigma;

ovrall = mean(D.acc(D.var == 1));
int = mean(D.acc(D.var == 2));
ord = mean(D.acc(D.var == 3));
temp = mean(D.acc(D.var == 4));

disp(['Overall decoding accuracy = ' num2str(ovrall)]);
disp(['Integrated decoding accuracy = ' num2str(int)]);
disp(['Order decoding accuracy = ' num2str(ord)]);
disp(['Timing decoding accuracy = ' num2str(temp)]);

%%Display figure with barplot
% figure
% barplot(D.var,D.acc,'split',D.var,'style_bold','leg',{'Overall' ,'Integrated','Order','Timing'},'leglocation','north');
% title(['# of iterations=',num2str(n),'  TEMP=',num2str(vtemp),'  ORD=',num2str(vord),'  INTER=',num2str(vinter),'  NOISE=',num2str(vnoise)]);
% %axis([0 7 0 1.2]);
% ylabel('z acc');
% xlabel('factor');

%average data across iterations for RSA
D.data = cat(3, D.data{:});
D.data = mean(D.data,3);

%to mean across runs before running dissimilarity
% dataMean(1,:) = mean(D.data(1:4:24,:));
% dataMean(2,:) = mean(D.data(2:4:24,:));
% dataMean(3,:) = mean(D.data(3:4:24,:));
% dataMean(4,:) = mean(D.data(4:4:24,:));
% dataDistances = squareform(pdist(squeeze(dataMean)));
%
% figure
% imagesc(dataDistances)
% colorbar

%to mean across distances calculated from each run
dataDistance(:,:,1) = squareform(pdist(squeeze(D.data(1:8,:)))); dataDistance(:,:,2) = squareform(pdist(squeeze(D.data(9:16,:))));
dataDistance(:,:,3) = squareform(pdist(squeeze(D.data(17:24,:)))); dataDistance(:,:,4) = squareform(pdist(squeeze(D.data(25:32,:))));
dataDistance(:,:,5) = squareform(pdist(squeeze(D.data(33:40,:)))); dataDistance(:,:,6) = squareform(pdist(squeeze(D.data(41:48,:))));

meanDataDistance = mean(dataDistance,3);
D.distance = meanDataDistance;

figure
imagesc(meanDataDistance)
colorbar

save([saveDir '\simulation_prepProd_' num2str(prepProd) '_temp_' num2str(vtemp) '_ord_' num2str(vord) '_int_' num2str(vinter) '_noise_' num2str(vnoise) '_n_' num2str(n) '.mat'], 'D')

%%% Multi-dimensional scaling
RDMs.RDM = meanDataDistance;
RDMs.name = 'simulation';
RDMs.color = [0 0 1];

userOptions = prepProdSimu_defineUserOptions;

% rsa.MDSConditions(RDMs, userOptions);


