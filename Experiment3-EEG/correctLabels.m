function [Data_Stim_VC, relChan_S1] = correctLabels (Data_Stim_VC, relChan_S1)

oldLabels = {'F7'	4; 'F3'	3; 'Fz'	2; 'F4'	30; 'F8' 31; 'FC1' 7; 'FC2'	29; 'T7' 9;...
    'C3' 8; 'Cz' 24; 'C4' 25; 'T8'	26; 'CP5' 11; 'CP6' 22; 'P7' 15; 'P3'	14; 'Pz' 13;...
    'P4' 19; 'P8' 20; 'O1' 16; 'O2'	18; 'TP9' 10; 'TP10' 21; 'Fp1'	1; 'Fp2' 32; 'FT9'	5;...
    'FT10'	27; 'FC5' 6; 'FC6'	28; 'CP1' 12; 'CP2'	23; 'Oz' 17};
oldPhyCh = [4 3 2 30 31 7 29 9 8 24 25 26 11 22 15 14 13 19 20 16 18 10 21 1 32 5 27 6 28 12 23 17]';
newLabels = {'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'T7';'C3';'Cz';'C4';'T8';'TP9';'CP5';'CP1';'CP2';'CP6';'TP10';'P7';'P3';'Pz';'P4';'P8';'PO9';'O1';'Oz';'O2';'PO10'};
oldTr = Data_Stim_VC.trial;

for tr =1:length(Data_Stim_VC.trial)
    for ch =1:size(Data_Stim_VC.trial{1,tr},1)
        Data_Stim_VC.trial{1,tr}(ch,:) = oldTr{1,tr}(oldPhyCh==ch,:);
        fprintf ('taking channel:')
        disp(find(oldPhyCh==ch));
    end
end
Data_Stim_VC.label = newLabels;

end