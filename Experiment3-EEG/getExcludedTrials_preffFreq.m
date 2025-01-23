function excludTr = getExcludedTrials_preffFreq

excludTr{1,2}=[2, 6, 8, 14, 18, 39, 57, 70, 85, 95, 97, 98, 113, 146, 172, 196];

excludTr{2,2}=[3, 4, 12, 16, 17, 18, 23, 29, 38, 57, 61, 73, 78, 89, 93, 101, 106, 108, 123, 128, 130, 135, 137, 158, 164, 169, 172, 173, 186, 195, 201];

excludTrT{3,1}=[80, 175];
excludTrT{3,2}=[35, 46, 84, 94, 102, 109, 112, 113, 125, 129, 130, 133, 140, 152, 168, 173, 180, 182, 187, 190, 191, 192, 194, 202, 205, 211, 212, 214, 215, 217];

allTrials =1:max([excludTrT{3,1}(end) excludTrT{3,2}(end)])+length(excludTrT{3,1});
allTrials2 = allTrials;
allTrials(excludTrT{3,1})=[];
allTrials(excludTrT{3,2})=[];
excludTr{3,2} = setdiff(allTrials2',allTrials','rows')';

excludTr{4,2}=[32, 64, 76, 105, 110, 111, 113, 117, 122, 123, 133, 134, 136, 145, 146, 151, 153, 158, 159, 165, 166, 180, 181, 187, 190, 191, 192, 197, 204, 221, 222, 223];

excludTr{5,2}=[22, 29, 34, 55, 56, 57, 64, 65, 67, 68, 72, 84, 89, 96, 98, 105, 107, 108, 109, 110, 116, 118, 125, 129, 131, 135, 137, 139, 140, 141, 145, 148, 149, 151, 154, 157, 159, 162, 163, 173, 174, 176, 177, 179, 183, 193, 196, 204, 209, 215, 218, 219, 220, 222, 224];

excludTr{6,2}=[5, 6, 7, 20, 23, 35, 39, 68, 88, 96, 109, 113, 119, 127, 137, 138, 141, 143, 144, 147, 151, 157, 163, 173, 184, 186, 197, 198, 200, 203, 204, 206, 208, 216, 220, 221];
               %Fp2 and O2 were interpolated for this participant
               
excludTrT{7,1}=[18, 29, 217];
excludTrT{7,2}=[4,55,207];
allTrials =1:max([excludTrT{7,1}(end) excludTrT{7,2}(end)])+length(excludTrT{7,1});
allTrials2 = allTrials;
allTrials(excludTrT{7,1})=[];
allTrials(excludTrT{7,2})=[];
excludTr{7,2} = setdiff(allTrials2',allTrials','rows')';

excludTrT{8,1}=24;
excludTrT{8,2}=[56, 84, 196, 200];
allTrials =1:max([excludTrT{8,1}(end) excludTrT{8,2}(end)])+length(excludTrT{8,1});
allTrials2 = allTrials;
allTrials(excludTrT{8,1})=[];
allTrials(excludTrT{8,2})=[];
excludTr{8,2} = setdiff(allTrials2',allTrials','rows')';

excludTrT{9,1}=[141,142];
excludTrT{9,2}=[24, 37, 116, 160, 187, 191, 193, 196];
allTrials =1:max([excludTrT{9,1}(end) excludTrT{9,2}(end)])+length(excludTrT{9,1});
allTrials2 = allTrials;
allTrials(excludTrT{9,1})=[];
allTrials(excludTrT{9,2})=[];
excludTr{9,2} = setdiff(allTrials2',allTrials','rows')';

excludTr{10,2}=[29, 43, 55, 61, 80, 81, 85, 86, 90, 98, 101, 110, 133, 142, 145, 150, 153, 163, 172, 187, 198, 200, 205, 208, 214, 223];
                %Cz was interpolated for this subject
                
excludTrT{11,1}=30;
excludTrT{11,2}=[7, 18, 26, 34, 37, 39, 49, 71, 77, 105, 108, 110, 112, 113, 114, 115, 116, 117, 122, 128, 131, 135, 139, 141, 144, 146, 147, 150, 157, 158, 165, 169, 178, 185, 191, 198, 199];
                 
allTrials =1:max([excludTrT{11,1}(end) excludTrT{11,2}(end)])+length(excludTrT{11,1});
allTrials2 = allTrials;
allTrials(excludTrT{11,1})=[];
allTrials(excludTrT{11,2})=[];
excludTr{11,2} = setdiff(allTrials2',allTrials','rows')';

excludTr{12,2}=[55, 57, 197];

end