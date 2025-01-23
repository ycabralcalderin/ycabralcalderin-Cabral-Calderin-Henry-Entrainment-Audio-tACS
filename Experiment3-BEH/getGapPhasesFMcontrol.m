function [act,phaseBin,carrier,targetsegments] = getGapPhasesFMcontrol(Stim)

act            = [];
phaseBin       = [];
carrier        = [];
targetsegments = [];
presPhase      = zeros(15,2);
phases         = 0:2*pi/15:2*pi-(2*pi/15);
sf             = 44100;

for tr = 1:length(Stim) %Loop through trials
    
    hil            = hilbert(Stim(tr).FreqModulator);
    phas           = wrapTo2Pi(angle(hil));
    GapStartSample = [];
    
    %Loop through Gaps
    %making sure that gap locations are correct!
    for s = 1:length(Stim(tr).Audio)-2
        if Stim(tr).Audio(s+1)==0 && Stim(tr).Audio(s)~=0
            GapStartSample = [GapStartSample s+1];
        end
    end
    
    %RECALCULATE SAMPLE AS MIDDLE
    gaplengthinsamples = (Stim(tr).stiminfo.gapDur/1000) * sf; % multiply by fs to convert from s to number of samples
    gapMidSample       = round(GapStartSample + (gaplengthinsamples/2));%[408864 520339 745739]

    %GapInfo(tr).GapPhase=phas(fix(Stim(tr).stiminfo.gaplocationsinsamples));
    GapInfo(tr).GapPhase = phas(gapMidSample);
  
    carrFreq             = getcarrierFreqbyGap(Stim(tr),phas(gapMidSample),gapMidSample); %GET CARRIER
    
    for g = 1:Stim(tr).stiminfo.numtargets
        for bin = 1:length(phases)-1
            if GapInfo(tr).GapPhase(g)>phases(bin)&&GapInfo(tr).GapPhase(g)<phases(bin+1)
                GapInfo(tr).GapPhaseBin(g) = bin;
            elseif GapInfo(tr).GapPhase(g)>phases(length(phases))
                GapInfo(tr).GapPhaseBin(g) = length(phases);
            end
        end
        presPhase(GapInfo(tr).GapPhaseBin(g),1) = presPhase(GapInfo(tr).GapPhaseBin(g))+1;
        %presPhase(GapInfo(tr).GapPhaseBin(g),2)=presPhase(GapInfo(tr).GapPhaseBin(g))+1;
    end
    currPhases = [GapInfo(tr).GapPhaseBin; Stim(tr).stiminfo.targetbins];
    act        = [act phas(gapMidSample)];
    carrier    = [carrier carrFreq];
    phaseBin   = [phaseBin currPhases];
    targetsegments = [targetsegments Stim(tr).stiminfo.targetsegments'];
end
end