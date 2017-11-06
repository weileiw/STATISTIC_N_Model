clc
clear all
close all

addpath('~/Dropbox/MOCM/WEILEI/myfunc')
load ~/Dropbox/MOCM/DATA/no3obs_90x180x24.mat
load ~/Dropbox/MOCM/DATA/transport_v4.mat
load NF2d_NN_4.mat
load M3d90x180x24v2.mat
load N2F_90x180.mat 
MSK = MSKS;
%no3obs = DIN;
ineg = find(NF2d(:)<=0);
NF2d(ineg) = 0.01;

%NO3 = no3obs(:,:,1);%.*MSK.PAC(:,:,1);
NO3 = NO3_insitu;
iocn = find(M3d(:));
tmp = M3d+nan;
tmp(iocn) = 0;
tmp = tmp(:,:,1);
NF2d_ave = tmp;

for ll = 1:length(index)
    tmp = cell2mat(index{ll});

    NF2d_ave(tmp(1),tmp(2)) ...
        = nanmean([NF2d(tmp(1),tmp(2)),...
        NF2d(tmp(1),tmp(2)+1),NF2d(tmp(1),tmp(2)-1),...  
        NF2d(tmp(1)+1,tmp(2)),NF2d(tmp(1)-1,tmp(2)),...
        NF2d(tmp(1)+1,tmp(2)+1),NF2d(tmp(1)+1,tmp(2)-1),...
        NF2d(tmp(1)-1,tmp(2)+1),NF2d(tmp(1)-1,tmp(2)-1)]);

end

NF2d = NF2d_ave;
keyboard
R = [0,0.1,0.2,0.6,1,2,20];
%R = log10space(log10(0.01),log10(10),8);
N2F_M = cell(length(R),1);
N2F_O = cell(length(R),1);

for ii = 1:length(R)-1
    for jj = 1:length(NO3(:))

        if ~isempty(D2d{jj})&NO3(jj)>R(ii)&NO3(jj)<=R(ii+1)

            N2F_M{ii} = cat(2,N2F_M{ii},num2cell(NF2d(jj)));
            N2F_O{ii} = cat(2,N2F_O{ii},D2d{jj});
        else
            continue
        end

    end
end

MOD = [];
LAB1 = [];
for ii = 1:length(R)-1

    %  plot(ii*ones(1,length(cell2mat(N2F_M{ii}))),log10(cell2mat(N2F_M{ii})),'ko');
    %  hold on
    tmp = log10(cell2mat(N2F_M{ii}));
    MM(ii) = mean(tmp(isfinite(tmp)));
    MM_se(ii) = std(tmp(isfinite(tmp)));%/sqrt(length(tmp(isfinite(tmp))));
    MM1(ii) = geomean(cell2mat(N2F_M{ii}));
    MM1_se(ii) = (std(cell2mat(N2F_M{ii})));
    MOD = [MOD,cell2mat(N2F_M{ii})];
    LAB1 = [LAB1,R(ii+1)*ones(1,length(cell2mat(N2F_M{ii})))];
end
%hold off
OBS = [];
LAB2 = [];

figure(1)
for ii = 1:length(R)-1

    %  plot(ii*ones(1,length(cell2mat(N2F_O{ii}))),log10(cell2mat(N2F_O{ii})),'ko');
    %  hold on
    tmp = log10((cell2mat(N2F_O{ii})));
    NN(ii) = mean(tmp(isfinite(tmp)));
    NN_se(ii) = std(tmp(isfinite(tmp)));%/sqrt(length(tmp(isfinite(tmp))));
    NN1(ii) = geomean((cell2mat(N2F_O{ii})));
    NN1_se(ii) = std((cell2mat(N2F_O{ii})));
    TT = cell2mat(N2F_O{ii});
    OBS = [OBS,cell2mat(N2F_O{ii})];
    LAB2 = [LAB2,R(ii+1)*ones(1,length(cell2mat(N2F_O{ii})))];
    %plot(ii*ones(1,length(TT)),log10(TT),'ko')
    %hold on
    %legend('off')
end
boxplot(log10(OBS*2),LAB2)
hold on
boxplot(log10(MOD),LAB1,'PlotStyle','compact')
hold on
L = MM_se;
U = MM_se;
e = errorbar([1:length(R)-1],MM,L,U,'-s',...
    'MarkerSize',5,'MarkerEdgeColor',...
    'red','MarkerFaceColor','red');
e.Color = 'red';
set(e,'linewidth',1)

%legend('Model')
%L = NN./NN_se;
%U = NN.*NN_se;
%e = errorbar([1:length(R)-1],log10(NN1),NN_se,'-d',...
%	     'MarkerSize',1,'MarkerEdgeColor',...
%	     'green','MarkerFaceColor','blue');
%e.Color = 'green';
%legend('Observation')
hold off


%figure(2)
%e = errorbar([1:length(R)-1],MM,MM_se,'s','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
%e.Color = 'red';
%hold on

xticks([1:length(R)])
xticklabels({'0-0.1','0.1-0.2','0.2-0.6','0.6-1','1-2','>2'});
ylim([-3,4])
yticks([-3:4]);
yticklabels({'0.001','0.01','0.1','1','10','100','1000','10000'});
xticklabel_rotate([],45,[])
legend('Model','Observation')
%text(1.1,2,'N = 189')
%text(2.1,2,'N = 123')
%text(3.1,2,'N = 167')
%text(4.1,2,'N = 30')
%text(5.1,2,'N = 80')
%text(5.5,1.5,'N = 30')

text(1.1,2,'N = 249')
text(2.1,2,'N = 111')
text(3.1,2,'N = 142')
text(4.1,2,'N = 63')
text(5.1,2,'N = 25')
text(5.5,1.5,'N = 33')

%e = errorbar([1:length(R)-1],log10(MM1),log10(MM1_se),'s','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
%e.Color = 'red';
%hold on
%e = errorbar([1:length(R)-1],log10(NN1),log10(NN1_se),'s','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
%e.Color = 'blue';
%hold off

xlabel('Observed NO_3 concentration (mmol m^-^3)')
ylabel('N_2 fixation (mmolN m^-^2 yr^-^1)')
