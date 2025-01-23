function carrFreq = getcarrierFreqbyGap(Stim,gapphase,gapsample)

%this function tries to calculate the carrier frequency of the audio stim
%at a given gap position

sf                 = 44100;
dur                = 20;
n                  = sf * dur;
t                  = (1:n) / sf;
carrier            = Stim.stiminfo.carrierfreq;
moddepth           = .67;
modfreq            = 2;
start_phases       = [0,pi/4,2*pi/4,3*pi/4,4*pi/4,5*pi/4,6*pi/4,7*pi/4]; %possible phases for the beginning of the stimulus
startphase         = start_phases(Stim.stiminfo.startphase);
%gapphase           = Stim.gapFinalmidPhase;%[1.92, 2.26, 3.66];
%gaplengthinsamples = round(Stim.gapDur * sf); % multiply by fs to convert from s to number of samples
%gapsample          = Stim.gaplocationsinsamples + (gaplengthinsamples/2);%[408864 520339 745739]
figure(5)
hold on
plot(t,carrier + (carrier*moddepth)*cos(2*pi*modfreq*t+startphase))
for ii = 1:length(gapphase)
    %     y = carrier + (carrier*moddepth)*cos(2*pi*modfreq*nearest(t,gaptime(ii)) + (startphase + gapphase(ii)));
    %     plot(gaptime(ii),y,'o')
    y        = carrier + (carrier*moddepth)*cos(2*pi*modfreq*t(gapsample(ii)) + (startphase));
    plot(t(gapsample(ii)),y,'o')
    carrFreq(ii) = y;
end
close figure 5
end