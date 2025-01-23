function  [preproc1_Stim] = preproc_EEG1_final_part3_prefFM(SubjFile,OutputPath,session)

diary(SubjFile)
diary on

load([OutputPath '_Data_Stim_VC_ICA.mat'],'Data_Stim_VC_ICA');
load([OutputPath '_Data_Stim_VC.mat'], 'Data_Stim_VC');

%LOAD ICA CLEANED DATA IF EXISTS
if exist([OutputPath '_Data_Stim_VC_ChanRep_ICAcleaned.mat'], 'file')
    load([OutputPath '_Data_Stim_VC_ChanRep_ICAcleaned.mat'], 'Data_Stim_VC_ChanRep_ICAcleaned')
else
    %perfect muscle correlation
    icaMuscle = corrPerfectMuscleArtf(Data_Stim_VC_ICA);
    
    %browse ICA results to identify the noise components
    cfg           = [];
    layout        = load('acticap-64ch-standard2.mat');
    cfg.layout    = layout.lay;
    %cfg.layout    = 'acticap-64ch-standard2.mat';
    cfg.viewmode  = 'component';
    cfg.compscale = 'local';
    ft_databrowser(cfg, Data_Stim_VC_ICA);
    
    %power analysis of ICA
    powComp (Data_Stim_VC_ICA);
    
    %look at the data and write down the components to be removed
    prompt   = {'Enter the IC to be removed from data';'';'';'';'';'';'';'';'';''};
    name     = 'ICA';
    numlines = 1;
    badComps = str2double(inputdlg(prompt,name,numlines));
    
    % remove noise components (using ICA data)
    cfg                             = [];
    cfg.component                   = badComps(badComps>=1); % to be removed component(s)
    Data_Stim_VC_ChanRep_ICAcleaned = ft_rejectcomponent(cfg, Data_Stim_VC_ICA, Data_Stim_VC);
    save([OutputPath '_Data_Stim_VC_ChanRep_ICAcleaned'], 'Data_Stim_VC_ChanRep_ICAcleaned', '-v7.3');
end

cfg                   = [];
cfg.lpfilter          = 'yes';                              % apply lowpass filter
cfg.lpfreq            = 30;                                  % lowpass at 30 Hz.
ICAcleanedLP          = ft_preprocessing(cfg,Data_Stim_VC_ChanRep_ICAcleaned);

cfg                                                  = [];
cfg.method                                           = 'summary';
%cfg.channel = 'EEG';
Data_LPStim_VC_ChanRep_ICAcleaned_visualCleaned = ft_rejectvisual (cfg, ICAcleanedLP);

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
    cfg.method        = 'distance';
    cfg.layout        = layout.lay;
    cfg.neighbourdist = 0.25;
    neighb            = ft_prepare_neighbours(cfg,ICAcleanedLP);
    badchannels       = ft_channelselection(Chan2Int,ICAcleanedLP.label);
    
    % interpolate bad channels
    cfg                         = [];
    cfg.badchannel              = badchannels;
    cfg.layout                  = layout.lay;
    cfg.neighbours              = neighb;
    ICAcleanedLP                = ft_channelrepair(cfg,ICAcleanedLP);
    
    cfg                                                  = [];
    cfg.method                                           = 'summary';
    %cfg.channel = 'EEG';
    Data_LPStim_VC_ChanRep_ICAcleaned_visualCleaned = ft_rejectvisual (cfg, ICAcleanedLP);
    
end

% check the data again
cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, Data_LPStim_VC_ChanRep_ICAcleaned_visualCleaned)

% %downsample
% cfg            = [];
% cfg.resamplefs = 100;  %%% define resample frequency
% Data_LPStim_VCds  = ft_resampledata(cfg,Data_LPStim_VC_ChanRep_ICAcleaned_visualCleaned);

%reref
cfg               = [];
cfg.channel       ='EEG';
cfg.implicitref   = 'FCz';
cfg.reref         = 'yes';
cfg.refchannel    = 'all';
cfg.refmethod     = 'avg';
Data_LPStim_VCrr = ft_preprocessing(cfg,Data_LPStim_VC_ChanRep_ICAcleaned_visualCleaned);

preproc1_Stim  = Data_LPStim_VCrr;
%save([OutputPath '_preprocessed_Stim_30HzLP_S' session],'preproc1_Stim','-v7.3');
end

function powComp (ICAdata)
% power spectrum component
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';%compute the power spectrum in all ICs
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:30;
freq = ft_freqanalysis(cfg, ICAdata);

nsubplots = length(ICAdata.label);
nbyn = round(sqrt(nsubplots));% sqrt(nsubplots) should not contain decimals,
type doc subplot
Nfigs = ceil(size(ICAdata.topo,2)/nsubplots);
tot = Nfigs*nsubplots;

rptvect = 1:size(ICAdata.topo,2);
rptvect = padarray(rptvect, [0 tot-size(ICAdata.topo,2)], 0,'post');
rptvect = reshape(rptvect,nsubplots,Nfigs)';
figure;
for r=1:size(rptvect,1)
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    k=0;
    for j=1:size(rptvect,2)
        
        if~(rptvect(r,j)==0)
            
            k=k+1;
            cfg=[];
            cfg.channel = freq.label(rptvect(r,j));
            %subplot(nbyn,nbyn,k);
            ft_singleplotER(cfg,freq);
        end
    end
end
end


function icaMuscle = corrPerfectMuscleArtf(ICAdata)
perfectMuscle = [-23.91,-25.43,-26.72,-27.50,-27.92,-28.11,-28.21,-28.20,-28.02,-27.64,-27.08,-26.49,-25.96,-25.40,-24.72,-23.95,-23.25,-22.73,-22.44,-22.35,-22.31,-22.21,-22.08,-21.99,-21.89,-21.78,-21.68,-21.60,-21.56,-21.52,-21.46,-21.40,-21.32,-21.20,-21.02,-20.89,-20.83,-20.77,-20.69,-20.66,-20.67,-20.64,-20.56,-20.43,-20.36,-20.37,-20.41,-20.43,-20.43,-20.44,-20.49,-20.55,-20.57,-20.58,-20.56,-20.54,-20.54,-20.55,-20.56,-20.58,-20.59,-20.58,-20.58,-20.61,-20.64,-20.63,-20.61,-20.61,-20.59,-20.55,-20.56,-20.57,-20.54,-20.53,-20.56,-20.59,-20.59,-20.60,-20.65,-20.70,-20.70,-20.66,-20.68,-20.74,-20.77,-20.77,-20.74,-20.71,-20.68,-20.67,-20.67,-20.68,-20.70,-20.71,-20.71,-20.69,-20.69,-20.71,-20.74,-20.73];

%data=ICAdata;
%downsample ICAdata
cfg            = [];
cfg.resamplefs = 500;  %%% define resample frequency
%cfg.detrend     = 'no';
data           = ft_resampledata(cfg,ICAdata);

Ntrials  = length(data.trial);
Nchans   = size(data.trial{1},1);
Nsamples = zeros(1,Ntrials);
for trial=1:Ntrials
    Nsamples(trial) = size(data.trial{trial},2);
end

dat = zeros(Nchans, sum(Nsamples));
for trial=1:Ntrials
    fprintf('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    dat(:,begsample:endsample) = data.trial{trial};
end

for k = 1:size(dat,1)
    [spec, f] = pwelch(dat(k,:),1024,0,1024,500);
    specSel = spec((f>1 &f<50));
    try
        corrR(k) = corr(perfectMuscle',log10(specSel));
    catch
        %this was here before
        %corrR(k) = corr(perfectMuscle(1:EEG.srate/500:end)',log10(specSel));
        error('modified this file for popout2 sampling rate of 1000, why did this not work with the other project?')
    end
    
    %     lowFreq(k) = mean(spectra(k,(freqs>1 & freqs <20)));
    %     highFreq(k) = mean(spectra(k,(freqs>20 & freqs <100)));
    lowFreq(k) = median(spec(f>1&f<20));
    highFreq(k)=  median(spec(f>20&f<80)); %before it was 100, but the data has been low-pass at 80Hz now
    
end
badSpect = highFreq./lowFreq > 2;
disp(corrR)
disp(badSpect)
% raw.spect = highFreq./lowFreq;
% raw.corr = corrR;
fprintf('\n Following Components have been marked :\n')
fprintf('%i ',find(corrR>0.6 | badSpect))
fprintf('\n Total of: %i \n',sum(corrR>0.6 | badSpect))
icaMuscle             = find(corrR>0.6 | badSpect); % find components that have a high correlation with "perfect" muscle components
ICAdata.icaMuscle     = icaMuscle;
ICAdata.icaMuscleCorr = corrR(corrR>0.6 | badSpect);
end
