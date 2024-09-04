function  [Data_Stim_VC_ICA] = preproc_EEG1_final_part2_prefFM(OutputPath)

% This function uses Fieldtrip to perform ICA on the previously
% pre-processed data

load([OutputPath '_Data_Stim_VC.mat'],'Data_Stim_VC');
load ([OutputPath 'relChan.mat'], 'relChan_S1');

% select only the EEG channels excluding the reference to run ICA in the
% % next step
% cfg              = [];
% cfg.channel      = {'EEG','-TP9', '-TP10'};
% dataMain  = ft_selectdata(cfg, visualCleaned);

% 7. - ICA for removing noise components
cfg                     = [];
cfg.method              = 'runica';  % this is the default and uses the implementation from EEGLAB
cfg.runica.extended     = 1;         % so the algorithm can also detect subgaussian sources of activity, such as line current and/or slow activity.
cfg.runica.stop         = 1e-7;      % to lower the criterion for stopping learning, thereby lengthening ICA training but possibly returning cleaner decompositions, particularly of high-density array data.
%cfg.runica.pca = 59 rank;
cfg.runica.pca          = relChan_S1; %;-1; % because we use average reference
Data_Stim_VC_ICA = ft_componentanalysis(cfg, Data_Stim_VC);
save([OutputPath '_Data_Stim_VC_ICA'],'Data_Stim_VC_ICA','-v7.3');

end

