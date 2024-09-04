%--EEG ANALYSIS USING RAW VALUES-----------------------------------------------

clear
close all
clc
%N = maxNumCompThreads
%LASTN = maxNumCompThreads(5) %don"t need

maxNumCompThreads(5);
% 0. Define some variables
mainDIR = uigetdir('path to raw data'); %pop up showing all folders to chose
eegDIR = fullfile(mainDIR,'EEG');
behDIR = fullfile(mainDIR,'BEH');
outDIR = fullfile(mainDIR,'ANA');

%addpath('/mnt/beegfs/users/yuranny.cabral/fieldtrip-20230503');
addpath('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/MATLAB_toolboxes/fieldtrip-20201128');
addpath('/mnt/beegfs/users/yuranny.cabral/FMtACSModDepth/EEGscripts');

% addpath('/hpc/users/yuranny.cabral/2019-0226-relentrai n/AnalysisScripts/');
% addpath('/hpc/users/yuranny.cabral/2019-0226-relentrain/AnalysisScripts/EEG1/');

%subjects = {'bestFreqEnt_24Feb2023';'bestFreqEnt_15March2023';'ECY06_1'}; %change subj names!
subjects = {'ATE26' 'BAO30' 'Bei24' 'CSA27' 'EGA23' 'Ei31' 'ESA20' 'LA03' 'lys16' 'MDS07' 'SUL03' 'WEC24'};

ft_defaults

% %creating variables for saving group data
% groupHitERP  = cell(numel(subjects),2);
% groupMissERP = cell(numel(subjects),2);

% %--------------LOOP THROUGH SUBJECTS
% for subj = 1:length(subjects)
% 
% %%PREPROCESSING
%     % try
%     SubjFile   = fullfile(eegDIR,[subjects{subj} '_EEG1']);
%     OutputPath = fullfile(outDIR,subjects{subj},subjects{subj});
%     fprintf('%s%s\n','... preprocessing file ', SubjFile)
%     %%    %close all force
%   %  if ~exist ([OutputPath '_preproc_Stim.mat'], 'file')
%         %1.PREPROCESSING before ICA
% %         [Data_Stim_VC, relChan_S1] = preproc_EEG1_final_part1_prefFM(SubjFile,OutputPath); %preproc before ICA
% %
% %         if strcmp(subjects{subj},'ESA20') || strcmp(subjects{subj},'lys16')
% %             %labels are messed up because the wrong workspace was used
% %             [Data_Stim_VC, relChan_S1] = correctLabels (Data_Stim_VC, relChan_S1);
% %         end
% %
% %         save([OutputPath 'relChan'],'relChan_S1');
% %         save([OutputPath '_Data_Stim_VC'],'Data_Stim_VC','-v7.3');
% % %
% %         %2. Preprocessing ICA
%    %    Data_LPStim_VC_ICA = preproc_EEG1_final_part2_prefFM(OutputPath);
% 
%         %3. Preprocessing after ICA
%         [preproc_Stim] = preproc_EEG1_final_part3_prefFM(SubjFile,OutputPath, '1');
%         save([OutputPath '_preproc_Stim'],'preproc_Stim','-v7.3');
% %     else
% %     end
% end

exclTinPrep = getExcludedTrials_preffFreq;

%FREQ ANA
for subj = 1:length(subjects)
    OutputPath = fullfile(outDIR,subjects{subj},subjects{subj});
    
    load([OutputPath '_preproc_Stim'],'preproc_Stim');
    
    %4. GET STIMULUS INFORMATION
    
    load(fullfile(mainDIR,'BEH',subjects{subj},'Stim',[subjects{subj} '_Stim_Main.mat']))
    exclT = exclTinPrep{subj,2};
    
    FMrateAll = nan(1,length(preproc_Stim.trial)+length(exclT));
    
    for s = 1:length(preproc_Stim.trial)+length(exclT)
        FMrateAll(s) = Stim(s).stiminfo.FM;
    end
    FMrateAll(exclT) = [];
    Stim(exclT) = [];   
       
    %6.Time freq analysis
    %run basic fft for virtual sensor data
    [RawMeanNorm, FMamp,phaseLockedData, fData] = runBasicFreqAnadiffFM_LAST(preproc_Stim,FMrateAll,Stim);
   % runBasicFreqAnadiffFM(preproc_Stim,FMrateAll);
    save([OutputPath '_basicFreqAna_VirtualChannel'],'RawMeanNorm', 'FMamp','phaseLockedData', 'fData','-v7.3');
     
     figHandles = findall(0,'Type','figure');
      
    % Loop through figures 2:end
    for i = 1:numel(figHandles)
        set(figHandles(i),'Units','inches');
        screenposition = get(figHandles(i),'Position');
        set(figHandles(i),...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print(figHandles(i), '-dpdf',[OutputPath 'IndivFigs.pdf','-append'])
        
    end
    
end
