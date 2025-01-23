function  [ratebyPhase, stimPhaseMod, meanRate, correctresp, RT, FalsA, HitbyCarrier, allresp] = get_hitRatebyPhase1500ms_1Cycle (Stim, responseMain, phaseBins)

% this is the case for EEG, for fMRI and tACS is different
fs                            = 44100;
correctresp                   = zeros(phaseBins,2);
RT(1:phaseBins)               = struct('RT',[]);
FalsA(1:length(responseMain)) = struct('Time',[]);
HitbyCarrier                  = zeros(length(responseMain),2);
allresp                       = [];

for trl = 1:length(responseMain) %loop across trials
     HitbyCarrier(trl,1)      = Stim(trl).stiminfo.carrierfreq;
     HitbyCarrier(trl,3)      = Stim(trl).stiminfo.numtargets;
    for target = 1:Stim(trl).stiminfo.numtargets %loop across targets
        targetbin                = Stim(trl).stiminfo.targetbins(target);
        correctresp(targetbin,1) = correctresp(targetbin,1)+1; %add a trial to the phase bin
        expectedT                = Stim(trl).stiminfo.gaplocationsinsamples(target)/fs;
        respTime                 = responseMain(trl).TimeKb-responseMain(trl).startStimGS;
        det                      = find(respTime>(expectedT+0.1) & respTime<(expectedT+1.5)); %find response in a 1.5 s window after the gap, add small time for reaction so to avoid too early responses???
        if det
            allresp                  = [allresp 1]; %1 if the target was detected
            correctresp(targetbin,2) = correctresp(targetbin,2)+1; %add a correct response to the phase bin
            RT(targetbin).RT         = [RT(targetbin).RT (respTime(det)-expectedT)];
            %CarF(targetbin).CF       = [CarF(targetbin).CF Stim(trl).Audio(fix(Stim(trl).stiminfo.gaplocationsinsamples(target)))];
            %NewStim2 = recStim(Stim(trl), fs);        
            %Get hit by carrier
            HitbyCarrier(trl,2)      =   HitbyCarrier(trl,2)+1;  
        else
            allresp                  = [allresp 0]; %0 if it was a miss
        end
    end
    %calculate false alarms with the opposite loop
    for kpress = 1:length(responseMain(trl).TimeKb) %loop across key presses
        expectedTall = Stim(trl).stiminfo.gaplocationsinsamples/fs;
        respTimeCurr = responseMain(trl).TimeKb(kpress)-responseMain(trl).startStimGS;
        if find(expectedTall<respTimeCurr & expectedTall>respTimeCurr-1.5) %if there was a gap in the preceding 1.5s window
        else
            if respTimeCurr>0 %to avoid taking responses that happened before the sound started
            FalsA(trl).Time        = [FalsA(trl).Time respTimeCurr]; %find response in a 1.5 s window after the gap, add small time for reaction so to avoid too early responses???
            end
        end
    end
end

disp(correctresp)
ratebyPhase = correctresp(:,2)./correctresp(:,1);
%ratebyPhase = [ratebyPhase; ratebyPhase]; % repeat for visualization

smoothrate   = smoothdata(ratebyPhase,'movmedian',2); % see paper from Molly NatComm
stimPhaseMod = max(smoothrate)-min(smoothrate); %stimulus-driven behavioral modulation
meanRate     = mean(ratebyPhase);
end