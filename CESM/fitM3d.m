load ~/Dropbox/MOCM/DATA/transport_v4.mat
%load ~/Dropbox/MOCM/DATA/omega4.mat
%[DAT,TXT] = xlsread('N2_fixation_clean.xlsx','No_Heter');
%[DAT,TXT] = xlsread('N2_fixation_clean.xlsx');
load Observed_N2F_NO3_90x180.mat
%lat = DAT(:,1);
%lon = DAT(:,2);
%ineg = find(lon<0);
%lon(ineg) = lon(ineg)+360;
%dat = DAT(:,3)*365/1000;
%no3_insitu = DAT(:,4);
grd = grid;
M3d = M3d;
D2d = cell(90,180);
D2d_no3_insitu = cell(90,180);

tt  = 1;
index{1} = num2cell([0,0]);

aa = length(dat);

for ii = 1:aa
    iflag = 0;
    Lat   = lat(ii);
    Lon   = lon(ii);

    for kk = 1:length(grd.xt)
        if Lon>=(grd.xt(kk)-1)&Lon<(grd.xt(kk)+1)
            a = kk;
        else
            continue
        end
    end

    for pp = 1:length(grd.yt)
        if Lat>=(grd.yt(pp)-1)&Lat<(grd.yt(pp)+1)
            b = pp;
        else
            continue
        end            
    end

    for rr = 1:length(index)
        temp = cell2mat(index{rr});
        if b == temp(1)&a == temp(2)
            D2d{b,a} = cat(2,D2d{b,a},num2cell(dat(ii)));
            D2d_no3_insitu{b,a} = cat(2,D2d_no3_insitu{b,a},...
                num2cell(no3_insitu(ii)));
            iflag = 1;
        else
            continue                
        end
    end

    if iflag == 0
        index{tt} = num2cell([b,a]);
        D2d{b,a} = cat(2,D2d{b,a},num2cell(dat(ii)));
        D2d_no3_insitu{b,a} = cat(2,D2d_no3_insitu{b,a},...
            num2cell(no3_insitu(ii)));
        tt = tt+1;
    end

end

iocn = find(M3d(:));
NO3_insitu = M3d+nan;
NO3_insitu(iocn) = 0;
NO3_insitu = NO3_insitu(:,:,1);

for ll = 1:length(index)

    tmp = cell2mat(index{ll});

    NO3_insitu(tmp(1),tmp(2)) ...
        = geomean(cell2mat(D2d_no3_insitu{tmp(1),tmp(2)}));
    STD(tmp(1),tmp(2)) ...
        = std(cell2mat(D2d_no3_insitu{tmp(1),tmp(2)}));

end

%DAT = DAT.*M3d; % at the very edge of land, the M3d mask 
% treats land as ocean.

fname = sprintf('N2F_90x180');
save(fname,'D2d','NO3_insitu','STD','index')




