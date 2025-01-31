%get tACS frequencies
clear
clc
close all

load('percentChange.mat')
load('/Users/yuranny.cabral/ownCloud - yuranny.cabral-calderin@ae.mpg.de@owncloud.gwdg.de/Paper_tACS_FMmodDepth/FinalScripts/OscillatorModel_amp','amp')
%model variable is strength x tacs/baseline x oscifreq x tACS freq

%first compute percent change to baseline
ampChange = squeeze((amp(:, 2,:,:) - amp(:, 1,:,:))./amp(:, 1,:,:));
oscPrefFreq      = 0.1:0.1:8;      %frequency of the ongoing oscillator in Hz

%assert(size(ampChange, 1) == length(percentChange), 'Mismatch between amp and percentChange dimensions!');


for subj =1:length(percentChange)

    %correlate tACS measured profile with model

    %correlate real to model
    for freq =1:size(ampChange,2)
        for strength =1:size(ampChange,1)
            SubjCorr(freq,strength) = corr(squeeze(ampChange(strength,freq,:)),percentChange(subj,:)');
        end
    end

    [ampChange2PLOT(subj),Pos] = max(SubjCorr,[],'all');
    [row,col] = ind2sub(size(SubjCorr),Pos);
    GroupRF(subj) = oscPrefFreq(row);
% %try getting the frequency by averaging across intensities
% %firt to fisher'z transformation
% SubjCorr_fishersZ   = atanh(SubjCorr);
% 
% 
%    [ampChange2PLOT(subj),Pos] = max(mean(SubjCorr_fishersZ,2));
%     [row,col] = ind2sub(size(SubjCorr_fishersZ),Pos);
%     GroupRF(subj) = oscPrefFreq(row);

    figure(1)
    subplot(2,12,subj)
    plot(zscore(squeeze(ampChange(col,row,:))))
    hold on
    plot(zscore(percentChange(subj,:))')
    title(['RF: ' num2str(GroupRF(subj))])

    figure(2)
    subplot(2,12,subj)
    histogram(SubjCorr,20)

    figure(3)
    subplot(2,12,subj)
    distributionPlot(SubjCorr','showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

    figure(4)
    subplot(2,12,subj)
    plot(oscPrefFreq,max(SubjCorr,[],2))
    hold on
    plot(GroupRF(subj), ampChange2PLOT(subj),'r*')
    xlabel ('Intrinsic oscillator frequency')
    ylabel ('r')
end

figure, distributionPlot(GroupRF','showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('predicted Resonance Frequency')
% xlabel('frequency')
% set(gca, 'XTick',[0.8 2 3.2 4.4])
% xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('Frequency')

%PLOT model-real data max correlation

figure, distributionPlot(ampChange2PLOT','showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('Correlation model to empirical data')
% xlabel('frequency')
% set(gca, 'XTick',[0.8 2 3.2 4.4])
% xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('Frequency')

%perform clustering
%cluster anaylsis
figure
eucD = pdist(GroupRF','euclidean');
clustTreeEuc = linkage(eucD,'average');
cophenet(clustTreeEuc,eucD)
[H,T] = dendrogram(clustTreeEuc,'Orientation','left','ColorThreshold','default');
set(H,'LineWidth',2)

figure
[~,T3] = dendrogram(clustTreeEuc,2,'Orientation','left','ColorThreshold','default');
clusterIDcommon = T3;
%PLOT

figure, subplot(2,2,1)
title(strcat('mean RF: ', num2str(mean(GroupRF(clusterIDcommon==1)))))
distributionPlot(percentChange(clusterIDcommon==1,:),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
hold on
plot(percentChange(clusterIDcommon==1,:)','--','color',[0.6 0.6 0.6],'LineWidth',1)
plot(mean(percentChange(clusterIDcommon==1,:))','g','LineWidth',2)
ylabel('Behavioral Entrainment')
ylim([-5 8])
subplot(2,2,3)
plotSpread(GroupRF(clusterIDcommon==1)')
ylim([0 5])

subplot(2,2,2)
title(strcat('mean RF: ', num2str(mean(GroupRF(clusterIDcommon==2)))))
distributionPlot(percentChange(clusterIDcommon==2,:),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
hold on
plot(percentChange(clusterIDcommon==2,:)','--','color',[0.6 0.6 0.6],'LineWidth',1)
plot(mean(percentChange(clusterIDcommon==2,:))','r','LineWidth',2)
ylabel('Behavioral Entrainment')
ylim([-5 8])
subplot(2,2,4)
plotSpread(GroupRF(clusterIDcommon==2)')
ylim([0 5])
% %Plotting figure with comparing across tasks
% %define subject and protocol info
% SUBJlist = {'ABI24' 'ATE26' 'BAO30' 'BBB19' 'BWC24' 'CSA27' 'EDI31' 'EGA23' 'EOE26' 'ESA20' 'HAE30' 'KWA10' 'LCS16' 'LPA03' 'MDS07' 'MUL19' 'OLDO03' 'PLU28' 'SFE21' 'SKN24' 'SRL03' 'TDN16' 'VEI09' 'ZJE04'};
% %EEGprefereFreList
% %%%% IN EEG PREFERED FREQ EXP:
% %ABI24=Bei24
% %EDI31=Ei31
% %LPA03=LA03
% %LCS16=lys16
% %SRL03=SUL03
% %BWC24=WEC24
% %subjects = {'ATE26' 'BAO30' 'Bei24' 'CSA27' 'EGA23' 'Ei31' 'ESA20' 'LA03' 'lys16' 'MDS07' 'SUL03' 'WEC24'};
% commonSubj = [1 2 3 5 6 7 8 10 13 14 15 21];
% clusterIDcommon = T3(commonSubj);

%plot behavior diffFM according to model
load('/Volumes/Projects/2021-0293-EntEcho/diffFMsFinal/Data/BEH/GroupData_diffFM_1500ms.mat')

%same as before but now for neural data
load('/Volumes/Projects/2021-0293-EntEcho/diffFMsFinal/diffFM_EEG_Group_Results.mat','GroupMeanNormV','GroupMeanNorm')

%correlate amplitude across conditions
% %define subject and protocol info
% SUBJlist = {'ABI24' 'ATE26' 'BAO30' 'BBB19' 'BWC24' 'CSA27' 'EDI31' 'EGA23' 'EOE26' 'ESA20' 'HAE30' 'KWA10' 'LCS16' 'LPA03' 'MDS07' 'MUL19' 'OLDO03' 'PLU28' 'SFE21' 'SKN24' 'SRL03' 'TDN16' 'VEI09' 'ZJE04'};
% %EEGprefereFreList
% %%%% IN EEG PREFERED FREQ EXP:
% %ABI24=Bei24
% %EDI31=Ei31
% %LPA03=LA03
% %LCS16=lys16
% %SRL03=SUL03
% %BWC24=WEC24
% %subjects = {'ATE26' 'BAO30' 'Bei24' 'CSA27' 'EGA23' 'Ei31' 'ESA20' 'LA03' 'lys16' 'MDS07' 'SUL03' 'WEC24'};
% commonSubj = [1 2 3 5 6 7 8 10 13 14 15 21];
reOrdertACS = [2 3 1 6 8 7 10 14 13 15 21 5];
[rho1,p1] = corr(groupfitamp(:,1),GroupMeanNormV(:,1,3),'type','Spearman')
[rho2,p2] = corr(groupfitamp(:,2),GroupMeanNormV(:,2,3),'type','Spearman')
[rho3,p3] = corr(groupfitamp(:,3),GroupMeanNormV(:,3,3),'type','Spearman')
[rho4,p4] = corr(groupfitamp(:,4),GroupMeanNormV(:,4,3),'type','Spearman')

[rho21,p21] = corr(groupfitamp(:,1),percentChange(reOrdertACS,1),'type','Spearman')
[rho22,p22] = corr(groupfitamp(:,2),percentChange(reOrdertACS,2),'type','Spearman')
[rho23,p23] = corr(groupfitamp(:,3),percentChange(reOrdertACS,3),'type','Spearman')
[rho24,p24] = corr(groupfitamp(:,4),percentChange(reOrdertACS,4),'type','Spearman')

[rho31,p31] = corr(GroupMeanNormV(:,1,3),percentChange(reOrdertACS,1),'type','Spearman')
[rho32,p32] = corr(GroupMeanNormV(:,2,3),percentChange(reOrdertACS,2),'type','Spearman')
[rho33,p33] = corr(GroupMeanNormV(:,3,3),percentChange(reOrdertACS,3),'type','Spearman')
[rho34,p34] = corr(GroupMeanNormV(:,4,3),percentChange(reOrdertACS,4),'type','Spearman')

[rho1_2,p1_2] = corr(groupfitamp(:,1),GroupMeanNorm(:,1,3),'type','Spearman')
[rho2_2,p2_2] = corr(groupfitamp(:,2),GroupMeanNorm(:,2,3),'type','Spearman')
[rho3_2,p3_2] = corr(groupfitamp(:,3),GroupMeanNorm(:,3,3),'type','Spearman')
[rho4_2,p4_2] = corr(groupfitamp(:,4),GroupMeanNorm(:,4,3),'type','Spearman')
[rho31_2,p31_2] = corr(GroupMeanNorm(:,1,3),percentChange(reOrdertACS,1),'type','Spearman')
[rho32_2,p32_2] = corr(GroupMeanNorm(:,2,3),percentChange(reOrdertACS,2),'type','Spearman')
[rho33_2,p33_2] = corr(GroupMeanNorm(:,3,3),percentChange(reOrdertACS,3),'type','Spearman')
[rho34_2,p34_2] = corr(GroupMeanNorm(:,4,3),percentChange(reOrdertACS,4),'type','Spearman')


%PLOT common subject entrainment amplitude by frequency

figure
for s =1:length(reOrdertACS)
    subplot(2,6,s)

    hold on
    plot(zscore(percentChange(reOrdertACS(s),:)))
    plot(zscore(groupfitamp(s,:)))
    plot(zscore(squeeze(GroupMeanNormV(s,:,3))))

end
