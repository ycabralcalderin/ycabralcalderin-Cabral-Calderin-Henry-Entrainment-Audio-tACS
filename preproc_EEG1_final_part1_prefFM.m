function  [Data_dmStim_VC, relChan] = preproc_EEG1_final_part1_prefFM(SubjFile,OutputPath)

%diary(SubjFile)
%diary on
if exist(OutputPath(1:end-5),'dir')
else
    mkdir([OutputPath(1:end-5) '/'])
end

%%%%%% Prep until before ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.-Defining trials
%%
cfg              = []; %cfg is a configuration structure
%%
cfg.trialfun     = 'trialfun_EEG1_diffFM';
cfg.headerfile   = [SubjFile '.vhdr'];
cfg.datafile     = [SubjFile '.eeg'];
cfg              = ft_definetrial(cfg);
trlstim_EEG      = cfg.trl;
save([OutputPath '_trlstim_EEG_S'],'trlstim_EEG');

%read the data (EEG data)
cfg              = [];
cfg.datafile     = [SubjFile '.eeg'];
cfg.channel      = 1:32;%{'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'T7';'C3';'Cz';'C4';'T8';'TP9';'CP5';'CP1';'CP2';'CP6';'TP10';'P7';'P3';'Pz';'P4';'P8';'PO9';'O1';'Oz';'O2';'PO10'};% {'EEG'};
cfg.continuous   = 'yes';
data_EEG         = ft_preprocessing(cfg);

cfg = [];
cfg.viewmode = 'butterfly' ;
ft_databrowser(cfg, data_EEG);

%2.-Cut data into trials
cfg           = [];
cfg.trl       = trlstim_EEG(trlstim_EEG (:, 4)== 1, :);
itis = diff(cfg.trl(:,2));
if itis(1)< 2 %overlapping smaples
    cfg.trl(1,:) = [];
end
Data_Stim   = ft_redefinetrial(cfg, data_EEG); %function that cuts the data

ft_databrowser(cfg, Data_Stim)

%3.Demean
cfg                 = [];
cfg.channel         = {'EEG'};
%cfg.continuous        = 'yes';
cfg.demean          = 'yes';
% cfg.lpfilter        = 'yes';                              % apply lowpass filter
% cfg.lpfreq          = 30;                                 % lowpass at 30 Hz.
Data_dmStim         = ft_preprocessing(cfg,Data_Stim);

%4.-reject trials visually
cfg                 = [];
cfg.method          = 'summary'; %different types of methods, eg channel
cfg.channel         = 'EEG';
Data_dmStim_VC      = ft_rejectvisual (cfg, Data_dmStim);%_RR); %use z-score to check what to remove.

%IN CASE there are BAD CHANNELS TO INTERPOLATE
interp = inputdlg('Y or N','Interp Electrodes?',1);

if strcmp(interp{1,1},'Y')
    
    prompt   = {'Enter the electrodes to interpolate in a raw';'';'';'';'';'';'';'';'';'';''};
    name     = 'Interpolate Electrodes';
    numlines = 1;
    Chan2Int = inputdlg(prompt,name,numlines);
    
    %interpolate bad channels, prepare neighbours
    layout            = load('acticap-64ch-standard2.mat');
    cfg               = [];
    cfg.method        = 'distance'; %define neighbours (0.2 cm close) for each electrode
    cfg.layout        = layout.lay;
    cfg.neighbourdist = 0.2;  %in cm
    neighb            = ft_prepare_neighbours(cfg,Data_dmStim_VC);
    badchannels       = ft_channelselection(Chan2Int,Data_dmStim_VC.label); %select noisy channels
    
    % interpolate bad channels
    cfg                    = [];
    cfg.badchannel         = badchannels;
    cfg.layout             = layout.lay;
    cfg.neighbours         = neighb;
    Data_dmStim_VC_ChanRep = ft_channelrepair(cfg,Data_dmStim_VC);
    Data_dmStim_VC         = Data_dmStim_VC_ChanRep;
    relChan                = length(Data_dmStim_VC.label)-length(cfg.badchannel); %number of independent channels for running ICA, this changes when interpolating channels
    
    %visual inspection
    cfg                 = [];  % use only default options
    cfg.channel         = 'EEG';
    ft_databrowser(cfg,Data_dmStim_VC_ChanRep);
    save([OutputPath '_Chan2Int'], 'Chan2Int')
else
    relChan             = length(Data_dmStim_VC.label); %number of independent channels for running ICA, this changes when interpolating channels
    
end

cfg           = [];
ft_databrowser(cfg, Data_dmStim_VC); %Check that the trial cutting makes sense

%5. detrend
cfg            = [];
%cfg.resamplefs = 500;  %%% define resample frequency
cfg.detrend   = 'yes';
%Data_LPStim_VC  = ft_resampledata(cfg,Data_LPStim_VC);
Data_dmStim_VC  = ft_preprocessing(cfg,Data_dmStim_VC);
end
