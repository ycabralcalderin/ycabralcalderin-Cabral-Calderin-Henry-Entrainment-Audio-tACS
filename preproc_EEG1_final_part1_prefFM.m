function  [hpData_LPStim_RR_VC, relChan] = preproc_EEG1_final_part1_prefFM(SubjFile,OutputPath)

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
cfg.continuous   = 'yes';
data_EEG         = ft_preprocessing(cfg);

% 2. -Highpass filter before anything???
cfg            = [];
cfg.channel    = {'EEG'};
cfg.continuous = 'yes';
cfg.hpfilter   = 'yes';                              % apply highpass filter
cfg.hpfreq     = 0.6;
%cfg.hpfilttype = 'fir';
%cfg.hpfiltord  = 5; %filter order 
%cfg.hpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
% cfg.hpfiltwintype = highpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
hpData_EEG     = ft_preprocessing(cfg,data_EEG);

cfg = [];
cfg.channel  = {'EEG'};
cfg.viewmode = 'butterfly' ;
ft_databrowser(cfg, data_EEG);
ft_databrowser(cfg, hpData_EEG)

% 3. -Cut the data according to the Stim triggers
%Get epochs from Stim track
% [StimTrack_EEG, trCurrStimTrack, trCurrGapTrack] = getStimTrackTriggers(data_EEG,trlstim_EEG);
% StimTrackTrials                                  = [trCurrStimTrack; trCurrGapTrack];
% 
% if trlstim_EEG(1,1)<0 %in case the EEG started late and we cannot use 1 s before stim, in this case the trial can be excluded later
%     %    trlstim_EEG(1,3)     = trlstim_EEG(1,3)-trlstim_EEG(1,1)-1;
%     %trlstim_EEG(1,1)     = 1;
%     %   trCurrStimTrack(1,1) = 1;
%     %  trCurrStimTrack(1,3) =trlstim_EEG(1,3);
%     trCurrStimTrack(1,:) = [];
%     trlstim_EEG(1,:)     = [];
%     disp('......WARNING NOT ENOUGH PRE-STIM TIME FOR TRIAL 1... TRIAL will be excluded')
% end

% save([OutputPath '_StimTrackDefTrials_S' session],'StimTrackTrials');
% save([OutputPath '_StimTrack_EEG_S' session], 'StimTrack_EEG');
% load([OutputPath '_StimTrackDefTrials_S' session],'StimTrackTrials');
% load([OutputPath '_StimTrack_EEG_S' session], 'StimTrack_EEG');

%Cut Data according to stimulus periods using the StimTrack Info
cfg           = [];
cfg.trl       = trlstim_EEG(trlstim_EEG (:, 4)== 1, :);
itis = diff(cfg.trl(:,2));
if itis(1)< 2 %overlapping smaples
    cfg.trl(1,:) = [];
end
hpData_Stim   = ft_redefinetrial(cfg, hpData_EEG); %function that cuts the data

ft_databrowser(cfg, hpData_Stim)
%DataStimt     = ft_redefinetrial(cfg, StimTrack_EEG);

%4.Filtering options options (high pass, low pass, and spectral dft filter
%for line noise
cfg                   = [];
cfg.channel           = {'EEG'};
cfg.continuous        = 'yes';
cfg.lpfilter          = 'yes';                              % apply lowpass filter
cfg.lpfreq            = 30;                                 % lowpass at 30 Hz.
% cfg.dftfilter         = 'yes';
% cfg.dftfreq           = 50; % to
% cfg.dftreplace        = 'neighbour';%'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation (default = 'zero')
% cfg.dftbandwidth      = 1;%bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz (default = [1 2 3])
% cfg.dftneighbourwidth = 1;%bandwidth of frequencies neighbouring line noise frequencies, applies to spectrum interpolation, in Hz (default = [2 2 2])
%cfg.padding           = 10;
hpData_LPStim         = ft_preprocessing(cfg,hpData_Stim);
% cfg = [];
% ft_databrowser(cfg, hpData_LPStim);

%StimTrackcuthPCut = ft_redefinetrial(cfg, StimTrackcuthP); %See if you need it. Not very necessary
% %-Baseline-correction options
% cfg                 = [];
% cfg.demean          = 'yes';
% %cfg.datafile     = [SubjFile '.eeg'];
% % cfg.baselinewindow  = [-0.2 0];
% dataDemean = ft_preprocessing(cfg,dataCut);

%inspect the data
%cfg = [];
%cfg.channel= 'EEG';
%ft_databrowser(cfg,hpData_LPStim);

% %5.-Re-referencing options       %can be done at any timepoint. also at the end
cfg               = [];
cfg.channel       ='EEG';
cfg.implicitref   = 'FCz';
cfg.reref         = 'yes';
%cfg.refchannel    = {'TP9' 'TP10'};
cfg.refchannel    = 'all';
cfg.refmethod     = 'avg';
hpData_LPStim_RR  = ft_preprocessing(cfg,hpData_LPStim);
% save([OutputPath '_hpData_LPStim_RRave_S' session],'hpData_LPStim_RR','-v7.3');
% %end

%VISUALIZE the data to get some impression of its quality
% cfg = [];  % use only default options
% %cfg.channel= 'EEG';
% ft_databrowser(cfg,hpData_LPStim);%_RR);

%6.-reject trials visually
cfg                 = [];
cfg.method          = 'trial'; %differnttyoes of methods, eg channel
cfg.channel         = 'EEG';
hpData_LPStim_RR_VC    = ft_rejectvisual (cfg, hpData_LPStim_RR);%_RR); %use z-score to check what to remove. 

cfg                 = [];
cfg.method          = 'summary'; %differnttyoes of methods, eg channel
cfg.channel         = 'EEG';
hpData_LPStim_RR_VC    = ft_rejectvisual (cfg, hpData_LPStim_RR);%_RR); %use z-score to check what to remove. 


% remove the trials that have artifacts from the trl
%cfg.trl([1],:) = [];

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
    neighb            = ft_prepare_neighbours(cfg,hpData_LPStim_RR_VC);
    badchannels       = ft_channelselection(Chan2Int,hpData_LPStim_RR_VC.label); %select noisy channels
    
    % interpolate bad channels
    cfg                         = [];
    cfg.badchannel              = badchannels;
    cfg.layout                  = layout.lay;
    cfg.neighbours              = neighb;
    hpData_LPStim_RR_VC_ChanRep = ft_channelrepair(cfg,hpData_LPStim_RR_VC);
    hpData_LPStim_RR_VC         = hpData_LPStim_RR_VC_ChanRep;
    relChan                     = length(hpData_LPStim_RR_VC.label)-length(cfg.badchannel); %number of independent channels for running ICA, this changes when interpolating channels
    
    %visual inspection
    cfg                 = [];  % use only default options
    cfg.channel         = 'EEG';
    ft_databrowser(cfg,hpData_LPStim_RR_VC_ChanRep);
save([OutputPath '_Chan2Int'], 'Chan2Int')
else
    relChan             = length(hpData_LPStim_RR_VC.label); %number of independent channels for running ICA, this changes when interpolating channels
    
end

cfg           = [];
ft_databrowser(cfg, hpData_LPStim_RR_VC); %Check that the trial cutting makes sense

% %Downsample
% cfg            = [];
% cfg.resamplefs = 500;  %%% define resample frequency
% %cfg.detrend   = 'no';
% hpData_LPStim_RR_VC  = ft_resampledata(cfg,hpData_LPStim_RR_VC);
end
