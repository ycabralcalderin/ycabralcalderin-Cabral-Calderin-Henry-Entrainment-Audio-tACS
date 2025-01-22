% Modelling Entrainment amplitude as a function of the missmatch between
% tACS and ongoing oscillator frequency Cabral-Calderin & Henry 2025, Plos
% Biology

% Original Script from Krause et al 2022 (fig4a.m), Plos Biology, https://doi.org/10.1371/journal.pbio.3001650
% Modified by Y. Cabral-Calderin

clc
clear
close all

START            = 0;  % seconds
STOP             = 2000; % seconds
PHASE_DIFFERENCE = 0; % degrees
pcts             = .1:.1:1;%[0.05, 0.30, 0.75]; %stim intensity %
scale_factor     = sqrt(0.2/1);
k                = pcts * scale_factor;
tACSfreqs        = [0.8 2 3.2 4.4]; %tACS frequencies in Hz
oscPrefFreq      = 0.1:0.1:8;      %frequency of the ongoing oscillator in Hz

%run the model to get the amplitude estiamtes
for f = 1:length(oscPrefFreq)

    freqMiss = tACSfreqs/oscPrefFreq(f); %this is the relf parameter, i.e., frequency relative to ongoing oscillation (multiple)

    for ftACS = 1:length(freqMiss)

        for ii = 1:length(k)

            [nostim_t, nostim_y] = tacs_YCCPlosBIOL(START, STOP, 0, 0, 1,oscPrefFreq(f));
            [tacs_t, tacs_y] = tacs_YCCPlosBIOL(START, STOP, k(ii), deg2rad(PHASE_DIFFERENCE), freqMiss(ftACS),oscPrefFreq(f));
            amp(ii, 2,f,ftACS) = rms(tacs_y(tacs_t > 40, 1)) * sqrt(2);
            amp(ii, 1,f,ftACS) = rms(nostim_y(nostim_t > 40, 1)) * sqrt(2);
        end
    end
end

save('/Users/yuranny.cabral/ownCloud - yuranny.cabral-calderin@ae.mpg.de@owncloud.gwdg.de/Paper_tACS_FMmodDepth/FinalScripts/OscillatorModel_amp','amp')

%run the model only for the values to plot in Fig4a
pcts             = [0.05, 0.30, 0.75]; %stim intensity %
k                = pcts * scale_factor;
oscPrefFreq      = [2 2.6 3.2];      %frequency of the ongoing oscillator in Hz
tACSfreqs        = [0.8 2 3.2 4.4]; %tACS frequencies in Hz
amp              = [];
PHASE_DIFFERENCE = 180; % degrees

for f = 1:length(oscPrefFreq)

    freqMiss = tACSfreqs/oscPrefFreq(f); %this is the relf parameter, i.e., frequency relative to ongoing oscillation (multiple)
    figure
    for ftACS = 1:length(freqMiss)

        for ii = 1:length(k)

            [nostim_t, nostim_y] = tacs_YCCPlosBIOL(START, STOP, 0, 0, 1,oscPrefFreq(f));
            [tacs_t, tacs_y] = tacs_YCCPlosBIOL(START, STOP, k(ii), deg2rad(PHASE_DIFFERENCE), freqMiss(ftACS),oscPrefFreq(f));
            amp(ii, 2,f,ftACS) = rms(tacs_y(tacs_t > 40, 1)) * sqrt(2);
            amp(ii, 1,f,ftACS) = rms(nostim_y(nostim_t > 40, 1)) * sqrt(2);
            ax = subplot(length(k), 4, ii*4-(4-ftACS));
            h = plot(ax, nostim_t, nostim_y(:,1), tacs_t, tacs_y(:, 1), 'LineWidth', 1);
            %    vline(10);
            xlim(ax, [0 40]);
            if ii == 1
                legend(ax, h, {'No-stim', 'tACS'});
            end
            title(sprintf('k = %0.0f%% Change: %0.02f%% tACSfreq: %0.01f oscFreq: %0.01f ', ...
                100*pcts(ii), ...
                100*(amp(ii, 2,f,ftACS) - amp(ii, 1,f,ftACS)) / amp(ii, 1,f,ftACS),...
                tACSfreqs(ftACS),...
                oscPrefFreq(f)));
        end
    end
end

%run the model only for the values to plot in Fig4b
pcts             = 0.1:.1:1; %stim intensity %
k                = pcts * scale_factor;
oscPrefFreq      = 2;      %frequency of the ongoing oscillator in Hz
tACSfreqs        = 1:.01:3;
amp              = [];
for f = 1:length(oscPrefFreq)

    freqMiss = tACSfreqs/oscPrefFreq(f); %this is the relf parameter, i.e., frequency relative to ongoing oscillation (multiple)
    for ftACS = 1:length(freqMiss)

        for ii = 1:length(k)

            [nostim_t, nostim_y] = tacs_YCCPlosBIOL(START, STOP, 0, 0, 1,oscPrefFreq(f));
            [tacs_t, tacs_y] = tacs_YCCPlosBIOL(START, STOP, k(ii), deg2rad(PHASE_DIFFERENCE), freqMiss(ftACS),oscPrefFreq(f));
            amp(ii, 2,f,ftACS) = rms(tacs_y(tacs_t > 40, 1)) * sqrt(2);
            amp(ii, 1,f,ftACS) = rms(nostim_y(nostim_t > 40, 1)) * sqrt(2);
        end
    end
end
percentChange = squeeze(100*(amp(:, 2,f,:) - amp(:, 1,f,:))./amp(:, 1,f,:));
colors2Plot = colormap('autumn');
figure
hold on
for c =1:length(colors2plot)
plot((tACSfreqs./oscPrefFreq-1), percentChange(c,:),'Color',colors2Plot(c*25-24,:),'LineStyle','-','Marker','o','LineWidth',2)
legend(num2str(pcts(c)))
end
ylabel('Percent Change')
xlabel('Frequency missmatch')
