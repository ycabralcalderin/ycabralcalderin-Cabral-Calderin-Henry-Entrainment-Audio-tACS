%--BEHAVIORAL ANALYSIS USING RAW VALUES-----------------------------------------------
%Last checked on 20231206
clear
close all
clc
% 0. Define some variables
mainDIR = uigetdir();
outDIR = fullfile(mainDIR,'ANA/');
addpath('/Users/yuranny.cabral/Documents/MATLAB/CircStat2012a')

%subjects = {'BAO30' 'CSA27' 'EGA23' 'EI31' 'LYS16' 'MDS07' 'SUL03' 'WEC24'};
subjects = {'ATE26' 'BAO30' 'BEI24' 'CSA27' 'EGA23' 'EI31' 'ESA20' 'LA03' 'LYS16' 'MDS07' 'SUL03' 'WEC24'};
% subjectNumberintACS = [2 3 1 6 8 7 10 14 13 15 21 5];
% clusterIDinmodel = T3(subjectNumberintACS);

%SUBJlisttACS = {'ABI24' 'ATE26' 'BAO30' 'BBB19' 'BWC24' 'CSA27' 'EDI31' 'EGA23' 'EOE26' 'ESA20' 'HAE30' 'KWA10' 'LCS16' 'LPA03' 'MDS07' 'MUL19' 'OLDO03' 'PLU28' 'SFE21' 'SKN24' 'SRL03' 'TDN16' 'VEI09' 'ZJE04'};
%EEGprefereFreList
%%%% IN EEG PREFERED FREQ EXP:
%ABI24=Bei24
%EDI31=Ei31
%LPA03=LA03
%LCS16=lys16
%SRL03=SUL03
%BWC24=WEC24



phaseBins = 15;
%preallocating some variables first
groupfinalTemp = nan(length(subjects),6); %to save best tempo information
groupfinalTemp2 = nan(length(subjects),6); %to save best tempo information
gapSize = nan(length(subjects),4);
%parameters for fitting the cosine function
X         = 0:2*pi/15:2*pi-2*pi/15; %this is my phase vector
X         = X+2*pi/15/2;
lag       = pi;
intercept = 0.5;
amp       = 0.5;
params    = [lag intercept amp];

%--------------- LOOP THROUGH SUBJECTS
for subj = 9:length(subjects)

    currSubj = subjects{subj};
    disp(currSubj)

    %%%% load raw data variables
    load([mainDIR '/' currSubj '/1/' currSubj '_responseMain_FMcontrol1.mat'],'responseMain');
    load([mainDIR '/' currSubj '/Stim/' currSubj '_Stim_Main.mat']);
    load([mainDIR '/' currSubj '/1/' currSubj '_spontaneousTempo1.mat']);
    load([mainDIR '/' currSubj '/1/' currSubj '_2ndspontaneousTempo1.mat']);

    groupfinalTemp(subj,1) = mean([responseTempo1(1).Tempo responseTempo1(2).Tempo responseTempo1(3).Tempo responseTempo2(1).Tempo responseTempo2(2).Tempo responseTempo2(3).Tempo]); %overall mean
    groupfinalTemp(subj,2) = std([responseTempo1(1).Tempo responseTempo1(2).Tempo responseTempo1(3).Tempo responseTempo2(1).Tempo responseTempo2(2).Tempo responseTempo2(3).Tempo]); %overall std
    groupfinalTemp(subj,3) = mean([responseTempo1(1).Tempo responseTempo1(2).Tempo responseTempo1(3).Tempo]); % mean1
    groupfinalTemp(subj,4) = std([responseTempo1(1).Tempo responseTempo1(2).Tempo responseTempo1(3).Tempo]); % std1
    groupfinalTemp(subj,5) = mean([responseTempo2(1).Tempo responseTempo2(2).Tempo responseTempo2(3).Tempo]); %mean2
    groupfinalTemp(subj,6) = std([responseTempo2(1).Tempo responseTempo2(2).Tempo responseTempo2(3).Tempo]); %std2

    %calculate preferred tempo based on median of ITI of last 25 taps
    groupfinalTemp2(subj,1)=mean([median(diff(responseTempo1(1).TimeKb )) median(diff(responseTempo1(2).TimeKb )) median(diff(responseTempo1(3).TimeKb ))...
        median(diff(responseTempo2(1).TimeKb )) median(diff(responseTempo2(2).TimeKb )) median(diff(responseTempo2(3).TimeKb ))]); %overall mean
    groupfinalTemp2(subj,2)=std([median(diff(responseTempo1(1).TimeKb )) median(diff(responseTempo1(2).TimeKb )) median(diff(responseTempo1(3).TimeKb ))...
        median(diff(responseTempo2(1).TimeKb )) median(diff(responseTempo2(2).TimeKb )) median(diff(responseTempo2(3).TimeKb ))]); %overall std
    groupfinalTemp2(subj,3)=mean([median(diff(responseTempo1(1).TimeKb )) median(diff(responseTempo1(2).TimeKb )) median(diff(responseTempo1(3).TimeKb ))]); % mean1
    groupfinalTemp2(subj,4)=std([median(diff(responseTempo1(1).TimeKb )) median(diff(responseTempo1(2).TimeKb )) median(diff(responseTempo1(3).TimeKb ))]); % std1
    groupfinalTemp2(subj,5)=mean([median(diff(responseTempo2(1).TimeKb )) median(diff(responseTempo2(2).TimeKb )) median(diff(responseTempo2(3).TimeKb ))]); %mean2
    groupfinalTemp2(subj,6)=std([median(diff(responseTempo2(1).TimeKb )) median(diff(responseTempo2(2).TimeKb )) median(diff(responseTempo2(3).TimeKb ))]); %std2


    % separate trials per FM condition
    FM = [0.8 2 3.2 4.4];
    idFM1 = [];
    idFM2 = [];
    idFM3 = [];
    idFM4 = [];

    for tr = 1:length(Stim)
        if Stim(tr).stiminfo.FM == FM(1)
            idFM1 = [idFM1 tr];
            if isnan(gapSize(subj,1))
                gapSize(subj,1)=Stim(tr).stiminfo.gapDur;
            end
        elseif Stim(tr).stiminfo.FM == FM(2)
            idFM2 = [idFM2 tr];
            if isnan(gapSize(subj,2))
                gapSize(subj,2)=Stim(tr).stiminfo.gapDur;
            end
        elseif Stim(tr).stiminfo.FM == FM(3)
            idFM3 = [idFM3 tr];
            if isnan(gapSize(subj,3))
                gapSize(subj,3)=Stim(tr).stiminfo.gapDur;
            end
        elseif Stim(tr).stiminfo.FM == FM(4)
            idFM4 = [idFM4 tr];
            if isnan(gapSize(subj,4))
                gapSize(subj,4)=Stim(tr).stiminfo.gapDur;
            end
        end
    end

    fmCode = [idFM1; idFM2; idFM3; idFM4];

    %loop trough FM conditions
    for currFM = 1:length(FM)    
        % 1. Calculate hit rate as a function of phase bin
        [groupData(subj).ratebyPhase(currFM,:)] = get_hitRatebyPhase1500ms_1Cycle(Stim(fmCode(currFM,:)), responseMain(fmCode(currFM,:)), phaseBins);
        [act,phaseBin,carrier,targetsegments] = getGapPhasesFMcontrol(Stim(fmCode(currFM,:)));

        if phaseBin(1,:)-phaseBin(2,:)==0 %double check that phases are correctly identified
        else
            error('Check wrong phases')
        end

        % 2. FIT cosine function to the behavioral modulation
        [groupData(subj).fit(currFM,:), groupData(subj).resnorm(currFM),groupData(subj).residual(currFM,:),groupData(subj).prefPhase(currFM) ] = fitCos2RatebyPhase_1cycle (groupData(subj).ratebyPhase(currFM,:)',X,params);
        yAll(currFM,:) = groupData(subj).fit(currFM,2) + groupData(subj).fit(currFM,3).*(cos(X + groupData(subj).fit(currFM,1))); % if you want to plot the predicted function
        disp(groupData(subj).prefPhase)
        groupData(subj).meanRate(currFM) = mean(groupData(subj).ratebyPhase(currFM,:));
    end

    groupData(subj).fit(currFM,:) = groupData(subj).fit(currFM,:)';

%     %DO some plotting
%     figure(subj)
%     XVal = (1:1:length(groupData(subj).ratebyPhase(currFM,:))); %this is just for plotting
%     subplot(2,2,1)
%     f    = 2;
%     %Amp  = 1;
%     ts   = 1/1000;
%     T    = .5;
%     t    = 0:ts:T;
%     yCos = cos(2*pi*f*t);
%     %make y positive
%     yCos = yCos+1;
%     yCos = yCos./max(yCos);
%     plot(phaseBins/length(yCos):phaseBins/length(yCos):phaseBins,yCos,'k')
%     iniColorAll = [0 0 1; 1 0 0; 0 1 0; 0 0 0];
%     hold on
%     for currFM = 1:4
%         iniColor=iniColorAll(currFM,:);
%         plot(XVal,groupData(subj).ratebyPhase(currFM,:),'.');
%         yNew = interp(yAll(currFM,:),5);
%         XValNew = interp(XVal,5);
%         hold on, plot(XValNew,yNew,'Color',iniColor,'LineStyle','-.')
%     end
%     phases = 0:2*pi/phaseBins:2*pi-2*pi/phaseBins;
%     meanphases = round(phases+2*pi/phaseBins/2,2);
%     xlabel ('Phase')
%     ylabel ('Hit rate')
%     
%     set(gca, 'XTick',1:1:phaseBins)
%     meanphases = meanphases;
%     xticklabels(meanphases)
%     
%     hold on
%     plot(1:1:phaseBins,yCos(1:phaseBins),'k')
%     title('Hit rate by FM Phase')
%     subplot(2,2,2)
%     circ_plot(groupData(subj).prefPhase','pretty')
%     subplot(2,2,3)
%     bar(groupData(subj).fit(:,2))
%     title('Hit rate')
%     xlabel('FM(Hz)')
%     ylabel('hit rate')
%     xticklabels([0.8 2 3.2 4.4])
%     subplot(2,2,4)
%     bar(groupData(subj).fit(:,3))
%     xticklabels([0.8 2 3.2 4.4])
%     title('sin fit amplitude')
%     ylabel('amplitude')
%     xlabel('FM(Hz)')
end
%Plotting Individual Data
figure
XVal = (1:1:length(groupData(subj).ratebyPhase(currFM,:))); %this is just for plotting
f    = 2;
%Amp  = 1;
ts   = 1/1000;
T    = .5;
t    = 0:ts:T;
yCos = cos(2*pi*f*t);
%make y positive
yCos = yCos+1;
yCos = yCos./max(yCos);

for subj=1:length(subjects)
    subplot(4,3,subj)
    plot(phaseBins/length(yCos):phaseBins/length(yCos):phaseBins,yCos,'k')
    iniColorAll = [230/255 0/255 148/255; 186/255 31/255 181/255; 102/255 0/255 161/255; 41/255 5/255 161/255];
    hold on
    for currFM = 1:4
        iniColor=iniColorAll(currFM,:);
        plot(XVal,groupData(subj).ratebyPhase(currFM,:),'.','Color',iniColor);
        yAll(currFM,:) = groupData(subj).fit(currFM,2) + groupData(subj).fit(currFM,3).*(cos(X + groupData(subj).fit(currFM,1))); % if you want to plot the predicted function
        yNew = interp(yAll(currFM,:),5);
        XValNew = interp(XVal,5);
        hold on, plot(XValNew,yNew,'Color',iniColor,'LineStyle','-.')
    end
    phases = 0:2*pi/phaseBins:2*pi-2*pi/phaseBins;
    meanphases = round(phases+2*pi/phaseBins/2,2);
    xlabel ('Phase')
    ylabel ('Hit rate')

    set(gca, 'XTick',1:1:phaseBins)
    meanphases = meanphases;
    xticklabels(meanphases)

    hold on
    plot(1:1:phaseBins,yCos(1:phaseBins),'k')
    title('Hit rate by FM Phase')
    ylim([0 1])
end

%Plot gap size per FM rata
figure
distributionPlot(gapSize,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
xticklabels([0.8 2 3.2 4.4])
title('Gap Size')
ylabel('Gap size')
xlabel('FM(Hz)')

for subj=1:length(groupData)
    Pphase(subj,:)=groupData(subj).prefPhase;
    Lagphase(subj,:)=groupData(subj).fit(:,1);
end

for subj=1:length(subjects)
    groupfitamp(subj,:)=groupData(subj).fit(:,3);
    [maxr(subj),prefR(subj)]=max(groupfitamp(subj,:));
    hR(subj,:)=groupData(subj).meanRate;
end

%exclude participants with mean HR <0.2 & >0.8, because best phase can not
%be properly estimated
keepSubj=find(mean(hR,2)>0.20&mean(hR,2)<0.9);

hT = FM(prefR)';
groupfinalTemp2F(:,1)=1./groupfinalTemp2(:,1); %convert SMT from s to Hz
[rho,normp]=corr(groupfinalTemp(:,1),hT,'type','Spearman');

figure 
subplot(2,4,1)
circ_plot(Pphase(keepSubj,1),'pretty')%,[],16,true,true,'linewidth',2,'color','r');
subplot(2,4,2)
circ_plot(Pphase(keepSubj,2),'pretty')%,[],16,true,true,'linewidth',2,'color','r');
subplot(2,4,5)
circ_plot(Pphase(keepSubj,3),'pretty')%,[],16,true,true,'linewidth',2,'color','r');
subplot(2,4,6)
circ_plot(Pphase(keepSubj,4),'pretty')%,[],16,true,true,'linewidth',2,'color','r');
title('Best Stimulus phase')
subplot(1,4,3)
distributionPlot(hR,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('Hit rates')
subplot(1,4,4)
distributionPlot(groupfitamp,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('modulation amplitude')

figure, distributionPlot(hT,'xValues',2,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)

%%RUN some Stats
%FIRST CHECK FOR NORMALITY
[normh(1,1),normp(1,1)]=kstest(hR(:,1));
[normh(1,2),normp(1,2)]=kstest(hR(:,2));
[normh(1,3),normp(1,3)]=kstest(hR(:,3));
[normh(1,4),normp(1,4)]=kstest(hR(:,4));

[normh(2,1),normp(2,1)]=kstest(groupfitamp(:,1));
[normh(2,2),normp(2,2)]=kstest(groupfitamp(:,2));
[normh(2,3),normp(2,3)]=kstest(groupfitamp(:,3));
[normh(2,4),normp(2,4)]=kstest(groupfitamp(:,4));

[normh(3,1),normp(3,1)]=kstest(gapSize(:,1));
[normh(3,2),normp(3,2)]=kstest(gapSize(:,2));
[normh(3,3),normp(3,3)]=kstest(gapSize(:,3));
[normh(3,4),normp(3,4)]=kstest(gapSize(:,4));


%VL
hRTable = table([hR(:,1); hR(:,2); hR(:,3); hR(:,4)], [ones(length(hR(:,1)),1)-.2; ones(length(hR(:,1)),1)+1; ones(length(hR(:,1)),1)+2.2; ones(length(hR(:,1)),1)+3.4],'VariableNames',...
    {'hR','Freq'});
hRlm = fitlm(hRTable,'hR~Freq');
figure, plot(hRlm)

groupfitampTable = table([groupfitamp(:,1); groupfitamp(:,2); groupfitamp(:,3); groupfitamp(:,4)], [ones(length(groupfitamp(:,1)),1)-.2; ones(length(groupfitamp(:,1)),1)+1; ones(length(groupfitamp(:,1)),1)+2.2; ones(length(groupfitamp(:,1)),1)+3.4],'VariableNames',...
    {'groupfitamp','Freq'});
groupfitamplm = fitlm(groupfitampTable,'groupfitamp~Freq');
figure, plot(groupfitamplm)

gapSizeTable = table([gapSize(:,1); gapSize(:,2); gapSize(:,3); gapSize(:,4)], [ones(length(gapSize(:,1)),1)-.2; ones(length(gapSize(:,1)),1)+1; ones(length(gapSize(:,1)),1)+2.2; ones(length(gapSize(:,1)),1)+3.4],'VariableNames',...
    {'gapSize','Freq'});
gapSizelm = fitlm(gapSizeTable,'gapSize~Freq');
figure, plot(gapSizelm)

%HR
%Omnibus test
[PhR,TABLEhR,STATShR] = friedman(hR);

%posthoc
[hRstats.P(1), hRstats.hRH(1), hRstats.s(1).hRSTATS] = signrank(hR(:,1),hR(:,2),'method','approximate');
[hRstats.P(2), hRstats.hRH(2), hRstats.s(2).hRSTATS] = signrank(hR(:,1),hR(:,3),'method','approximate');
[hRstats.P(3), hRstats.hRH(3), hRstats.s(3).hRSTATS] = signrank(hR(:,1),hR(:,4),'method','approximate');
[hRstats.P(4), hRstats.hRH(4), hRstats.s(4).hRSTATS] = signrank(hR(:,2),hR(:,3),'method','approximate');
[hRstats.P(5), hRstats.hRH(5), hRstats.s(5).hRSTATS] = signrank(hR(:,2),hR(:,4),'method','approximate');
[hRstats.P(6), hRstats.hRH(6), hRstats.s(6).hRSTATS] = signrank(hR(:,3),hR(:,4),'method','approximate');
%compute bonferroni holms
for pit=1:length(hRstats.P)
hRstats.bonfHolms(pit)=hRstats.P(pit).*(sum(hRstats.P>hRstats.P(pit))+1);
end
hRstats.bonfHolms = tmpP(tmpP(:,1),3)';
hRstats.bonfP = hRstats.P*length(hRstats.P);

%AMP
%Omnibus test
[Pamp,TABLEamp,STATSamp] = friedman(groupfitamp);

%posthoc
[groupfitampstats.P(1), groupfitampstats.groupfitampH(1), groupfitampstats.s(1).groupfitampSTATS] = signrank(groupfitamp(:,1),groupfitamp(:,2),'method','approximate');
[groupfitampstats.P(2), groupfitampstats.groupfitampH(2), groupfitampstats.s(2).groupfitampSTATS] = signrank(groupfitamp(:,1),groupfitamp(:,3),'method','approximate');
[groupfitampstats.P(3), groupfitampstats.groupfitampH(3), groupfitampstats.s(3).groupfitampSTATS] = signrank(groupfitamp(:,1),groupfitamp(:,4),'method','approximate');
[groupfitampstats.P(4), groupfitampstats.groupfitampH(4), groupfitampstats.s(4).groupfitampSTATS] = signrank(groupfitamp(:,2),groupfitamp(:,3),'method','approximate');
[groupfitampstats.P(5), groupfitampstats.groupfitampH(5), groupfitampstats.s(5).groupfitampSTATS] = signrank(groupfitamp(:,2),groupfitamp(:,4),'method','approximate');
[groupfitampstats.P(6), groupfitampstats.groupfitampH(6), groupfitampstats.s(6).groupfitampSTATS] = signrank(groupfitamp(:,3),groupfitamp(:,4),'method','approximate');
%compute bonferroni holms
for pit=1:length(groupfitampstats.P)
groupfitampstats.bonfHolms(pit)=groupfitampstats.P(pit).*(sum(groupfitampstats.P>groupfitampstats.P(pit))+1);
end
groupfitampstats.bonfP = groupfitampstats.P*length(groupfitampstats.P);

%GAPSIZE
[PgapSize,TABLEgapSize,STATSgapSize] = friedman(gapSize);

%posthoc
[gapSizestats.P(1), gapSizestats.gapSizeH(1), gapSizestats.s(1).gapSizeSTATS] = signrank(gapSize(:,1),gapSize(:,2),'method','approximate');
[gapSizestats.P(2), gapSizestats.gapSizeH(2), gapSizestats.s(2).gapSizeSTATS] = signrank(gapSize(:,1),gapSize(:,3),'method','approximate');
[gapSizestats.P(3), gapSizestats.gapSizeH(3), gapSizestats.s(3).gapSizeSTATS] = signrank(gapSize(:,1),gapSize(:,4),'method','approximate');
[gapSizestats.P(4), gapSizestats.gapSizeH(4), gapSizestats.s(4).gapSizeSTATS] = signrank(gapSize(:,2),gapSize(:,3),'method','approximate');
[gapSizestats.P(5), gapSizestats.gapSizeH(5), gapSizestats.s(5).gapSizeSTATS] = signrank(gapSize(:,2),gapSize(:,4),'method','approximate');
[gapSizestats.P(6), gapSizestats.gapSizeH(6), gapSizestats.s(6).gapSizeSTATS] = signrank(gapSize(:,3),gapSize(:,4),'method','approximate');
for pit=1:length(gapSizestats.P)
gapSizestats.bonfHolms(pit)=gapSizestats.P(pit).*(sum(gapSizestats.P>gapSizestats.P(pit))+1);
end
gapSizestats.bonfP = gapSizestats.P*length(gapSizestats.P);

[p1,z1]=circ_rtest(PphaseNew(:,1));
[p2,z2]=circ_rtest(PphaseNew(:,2));
[p3,z3]=circ_rtest(PphaseNew(:,3));
[p4,z4]=circ_rtest(PphaseNew(:,4));

%load (fullfile(mainDIR,'GroupData_diffFM_1500ms'))
%Parametric Hotelling paired sample test for equal angular means
for var1 =1:4
    for var2 = 1:4
 [pval(var1,var2), F(var1,var2)] = circ_htest(PphaseNew(:,var1),PphaseNew(:,var2));
    end
end

%correct for multiple comparisons
Bonfpval = pval.*6;

%[pval, table] = circ_wwtest([PphaseNew(:,1);PphaseNew(:,2);PphaseNew(:,3);PphaseNew(:,4)],[ones(length(PphaseNew(:,1)),1);ones(length(PphaseNew(:,1)),1)+1;ones(length(PphaseNew(:,1)),1)+2;ones(length(PphaseNew(:,1)),1)+3])

%plot separated by cluster according to model

save (fullfile(mainDIR,'GroupData_diffFM_1500ms'))