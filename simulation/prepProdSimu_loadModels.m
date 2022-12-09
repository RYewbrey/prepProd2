function models = prepProdSimu_loadModels(modelDir)

%load representational models
cd(modelDir)
files = dir('*.mat'); %loads all .mat files in modelDir
for i=1:length(files) %make sure modelDir only contains model files
    load(files(i).name)
    models.(D.modelName) = D;
end%for files
clear('D')