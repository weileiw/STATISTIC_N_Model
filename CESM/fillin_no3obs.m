addpath('~/Dropbox/MOCM/WEILEI/myfunc')
load ~/Dropbox/MOCM/DATA/transport_v4.mat
load ~/Dropbox/MOCM/DATA/no3obs_90x180x24.mat
[DAT,TXT]  = xlsread('N2_fixation_clean.xlsx');
lat        = DAT(:,1);
lon        = DAT(:,2);
no3obs     = no3obs(:,:,1);
ineg       = find(lon<0);
lon(ineg)  = lon(ineg)+360;
dat        = DAT(:,3)*365/1000;
izero      = find(dat(:)==0);
dat(izero) = 0.001;
no3_insitu = DAT(:,4);
izero      = find(no3_insitu(:)==0);
no3_insitu(izero) = 0.001;
grd = grid;

%[X,Y,Z] = meshgrid(1:2:359,-89:2:89,grd.zt);
%[mu,var,n,Q] = bin3d(lon,lat,18*ones(length(lon),1),...
    %dat,(dat*0.5),X,Y,Z);
%NF_obs = M3d+nan;
%iwet = find(M3d(:));
%NF_obs(iwet) = 0;
%ifix = find(mu(:)>0);

%NF_obs(ifix) = mu(ifix);
%NF_obs = NF_obs(:,:,1);
%fname = sprintf('NF_obs');
%save(fname,'NF_obs');

%keyboard;
for kk = 1:length(dat)
    for ii = 1:length(grd.xt)
        for jj = 1:length(grd.yt)

            Lat = lat(kk);
            Lon = lon(kk);

            if (Lon>=(grd.xt(ii)-1)&Lon<(grd.xt(ii)+1)&...
                    Lat>=(grd.yt(jj)-1)&Lat<(grd.yt(jj)+1)&...
                    isnan(no3_insitu(kk)))
                no3_insitu(kk) = no3obs(jj,ii);
            else
                continue
            end
        end
    end
end

inan = find(isnan(no3_insitu));
lat(inan) = [];
lon(inan) = [];
dat(inan) = [];
no3_insitu(inan) = [];

%fname = sprintf('Observed_N2F_NO3_90x180');
%save(fname,'lat','lon','dat','no3_insitu')
