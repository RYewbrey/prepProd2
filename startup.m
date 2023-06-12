if strncmp(computer,'PC',2) %windows PC (at home or office)
    
    addpath(genpath('Z:/rhys/prepProd2/matlab')); %Adjust! loaded with subdirectories (genpath command)
    addpath('Z:/toolboxes/spm12');
    addpath(genpath('Z:/toolboxes/tools')); %joern's extensions for spm
    addpath(genpath('Z:/toolboxes/userfun')); %joern's util tools (open source)
    addpath(genpath('Z:/toolboxes/region-master')); %joern's region toolbox for spm
    addpath(genpath('Z:/toolboxes/spm12/toolbox/suit')); %SUIT Cerebellum
    addpath(genpath('Z:/toolboxes/spm12/toolbox/DARTEL')); %DARTEL deformation for suit reslice
    addpath(genpath('Z:/toolboxes/spm12/toolbox/OldSeg')); %for suit reslice
    addpath(genpath('Z:/toolboxes/permutest')); %permutest for crossSection analysis
    addpath('Z:/toolboxes/rsatoolbox_matlab'); %RSA toolbox
    addpath(genpath('Z:/toolboxes/pcm_toolbox')); %PCM toolbox
    
else  %blueBear computing resources
    
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/rhys/prepProd2/matlab')); %Adjust! loaded with subdirectories (genpath command)
    addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12');
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/tools')); %joern's extensions for spm
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/userfun')); %joern's util tools (open source)
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/region-master')); %joern's region toolbox for spm
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/suit')); %SUIT cerebellum toolbox
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/DARTEL')); %DARTEL transformation
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/OldSeg')); %For SUIT reslicing
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/permutest')); %permutest for crossSection analysis
    addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/rsatoolbox_matlab'); %RSA toolbox
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/pcm_toolbox')); %PCM toolbox
    addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/surfing')); %surfing, used in RSA toolbox
    
end
