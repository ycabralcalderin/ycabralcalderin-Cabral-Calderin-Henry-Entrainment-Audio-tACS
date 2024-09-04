%from EEG results
%1 row: EEG preferred rate from vector length
%2 row: EEG preferred rate from amplitude
%3 BEH to FM preferred rate
%4 tACS best Freq
%5 computational model
Group_prefRate = [3 4 2 1 3 4 4 3 3 2 4 4;
                  3 4 3 4 2 4 4 2 3 4 4 4;
                  4 3 3 3 3 2 2 2 4 4 4 3;
                  1 1 2 1 1 1 1 2 3 2 4 2]';

figure, plot (Group_prefRate)

[rho1,p1]=corr(Group_prefRate(:,1),Group_prefRate(:,2),'type','Spearman')
[rho2,p2]=corr(Group_prefRate(:,1),Group_prefRate(:,3),'type','Spearman')
[rho3,p3]=corr(Group_prefRate(:,1),Group_prefRate(:,4),'type','Spearman')

[rho4,p4]=corr(Group_prefRate(:,2),Group_prefRate(:,3),'type','Spearman')
[rho5,p5]=corr(Group_prefRate(:,2),Group_prefRate(:,4),'type','Spearman')
[rho6,p6]=corr(Group_prefRate(:,3),Group_prefRate(:,4),'type','Spearman')