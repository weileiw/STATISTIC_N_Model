clc 
clear all
close all

load NF2d_OPT_N2P_tanh.mat
load ~/DATA/tempobs_90x180x24.mat 
load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat

TEMP = tempobs(:,:,1);

%iposi = and(NF2d(:)>0,TEMP(:)>0); % use only positive numbers.
iposi = find(NF2d(:)>0);
NF  = log10(NF2d(iposi));
TEMP = TEMP(iposi);

NF(imag(NF)~=0) = nan;    % change imaginary number to nan.
TEMP(imag(TEMP)~=0) = nan;  % change imaginary number to nan.
% get rid of nan and infinity numbers.
ivalue = and(~or(isnan(NF),isnan(TEMP)),~or(isinf(NF), isinf(TEMP)));

NF  = NF(ivalue);
TEMP = TEMP(ivalue);

NF_z  = (NF-mean(NF))/std(NF);
TEMP_z = (TEMP-mean(TEMP))/std(TEMP);

ATL_srf = MSKS.ATL(:,:,1);
ATL_srf = ATL_srf(iposi);
PAC_srf = MSKS.PAC(:,:,1);
PAC_srf = PAC_srf(iposi);
IND_srf = MSKS.IND(:,:,1);
IND_srf = IND_srf(iposi);
MED_srf = MSKS.MED(:,:,1);
MED_srf = MED_srf(iposi);
ARC_srf = MSKS.ARC(:,:,1);
ARC_srf = ARC_srf(iposi);

NF  = NF_z;
TEMP = TEMP_z;

% apply masks to fixation.
NF_ATL = NF.*ATL_srf(ivalue);
TEMP_ATL = TEMP.*ATL_srf(ivalue);

NF_PAC = NF.*PAC_srf(ivalue);
TEMP_PAC = TEMP.*PAC_srf(ivalue);

NF_IND = NF.*IND_srf(ivalue);
TEMP_IND = TEMP.*IND_srf(ivalue);

NF_ARC = NF.*ARC_srf(ivalue);
TEMP_ARC = TEMP.*ARC_srf(ivalue);

% ==========================================================================
fprintf('global ocean regression \n')
% linear regression 
RG_GLOB = fitlm(TEMP,NF,'linear');
fprintf('R^2 = %0.3f; \n',  RG_GLOB.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_GLOB.Coefficients.pValue(2))
% linear regression
[p,ErrorEst] = polyfit(TEMP,NF,1);
[NF_fit,delta] = polyval(p,TEMP,ErrorEst);

figure()
plot(TEMP,NF,'+',...
  TEMP,NF_fit,'g*',...
  TEMP,NF_fit+2*delta,'r<',...
  TEMP,NF_fit-2*delta,'r<')
title('linear regression for global ocean')

% quadratic regression
[p,ErrorEst] = polyfit(TEMP,NF,2);
[NF_fit,delta] = polyval(p,TEMP,ErrorEst);
[r2,rms] = rsquare(NF,NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(TEMP,NF,'+',...
  TEMP,NF_fit,'g*',...
  TEMP,NF_fit+2*delta,'r<',...
  TEMP,NF_fit-2*delta,'r<')
title('Quadratic regression for global ocean')

% ========================================================================
fprintf('Atlantic ocean regression \n')
% select data only in the Atlantic ocean for regression. 
ido = find(NF_ATL~=0);
RG_ATL = fitlm(TEMP_ATL(ido),NF_ATL(ido));

fprintf('R^2 = %0.3f; \n',  RG_ATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_ATL.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(TEMP_ATL(ido),NF_ATL(ido),2);
[NF_fit,delta] = polyval(p,TEMP_ATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_ATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(TEMP_ATL(ido),NF_ATL(ido),'+',...
  TEMP_ATL(ido),NF_fit,'g*',...
  TEMP_ATL(ido),NF_fit+2*delta,'r<',...
  TEMP_ATL(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Atlantic ocean')

%==========================================================================
fprintf('Pacific ocean regression \n')
% select data only in the Pacific ocean for regression. 
ido = find(NF_PAC~=0);
RG_PAC = fitlm(TEMP_PAC(ido),NF_PAC(ido));

fprintf('R^2 = %0.3f; \n',  RG_PAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_PAC.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(TEMP_PAC(ido),NF_PAC(ido),2);
[NF_fit,delta] = polyval(p,TEMP_PAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_PAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(TEMP_PAC(ido),NF_PAC(ido),'+',...
  TEMP_PAC(ido),NF_fit,'g*',...
  TEMP_PAC(ido),NF_fit+2*delta,'r<',...
  TEMP_PAC(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Pacific ocean')

%===========================================================================
fprintf('Indian ocean regression \n')
% select data only in the Indian ocean for regression. 
ido = find(NF_IND~=0);
RG_IND = fitlm(TEMP_IND(ido),NF_IND(ido));

fprintf('R^2 = %0.3f; \n',  RG_IND.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_IND.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(TEMP_IND(ido),NF_IND(ido),2);
[NF_fit,delta] = polyval(p,TEMP_IND(ido),ErrorEst);
[r2,rms] = rsquare(NF_IND(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(TEMP_IND(ido),NF_IND(ido),'+',...
  TEMP_IND(ido),NF_fit,'g*',...
  TEMP_IND(ido),NF_fit+2*delta,'r<',...
  TEMP_IND(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Indian ocean')

%===========================================================================
fprintf('Arctic ocean regression \n')
% select data only in the Arctic ocean for regression. 
ido = find(NF_ARC~=0);
RG_ARC = fitlm(TEMP_ARC(ido),NF_ARC(ido));

fprintf('R^2 = %0.3f; \n',  RG_ATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n \n \n',RG_ATL.Coefficients.pValue(2))

% ============================================================================
% create masks for north and south ocean sections.
NATL = MSKS.ATL;
NATL(1:45,:,:) = 0;
SATL = MSKS.ATL;
SATL(46:end,:,:) = 0;

NPAC = MSKS.PAC;
NPAC(1:45,:,:) = 0;
SPAC = MSKS.PAC;
SPAC(46:end,:,:) = 0;

NATL_srf = NATL(:,:,1);
NATL_srf = NATL_srf(iposi);
SATL_srf = SATL(:,:,1);
SATL_srf = SATL_srf(iposi);
NPAC_srf = NPAC(:,:,1);
NPAC_srf = NPAC_srf(iposi);
SPAC_srf = SPAC(:,:,1);
SPAC_srf = SPAC_srf(iposi);

% apply masks to fixation.
NF_NATL = NF.*NATL_srf(ivalue);
TEMP_NATL = TEMP.*NATL_srf(ivalue);
NF_SATL = NF.*SATL_srf(ivalue);
TEMP_SATL = TEMP.*SATL_srf(ivalue);

NF_NPAC = NF.*NPAC_srf(ivalue);
TEMP_NPAC = TEMP.*NPAC_srf(ivalue);
NF_SPAC = NF.*SPAC_srf(ivalue);
TEMP_SPAC = TEMP.*SPAC_srf(ivalue);

% =====================================================================
fprintf('N. Atlantic ocean regression \n')
% select data only in the N. Atlantic ocean for regression. 
ido = find(NF_NATL~=0);
RG_NATL = fitlm(TEMP_NATL(ido),NF_NATL(ido));
fprintf('R^2 = %0.3f; \n',  RG_NATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NATL.Coefficients.pValue(2))

[p,ErrorEst] = polyfit(TEMP_NATL(ido),NF_NATL(ido),1);
[NF_fit,delta] = polyval(p,TEMP_NATL(ido),ErrorEst);

figure()
plot(TEMP_NATL(ido),NF_NATL(ido),'k+',...
  TEMP_NATL(ido),NF_fit,'r*',...
  TEMP_NATL(ido),NF_fit+2*delta,'y<',...
  TEMP_NATL(ido),NF_fit-2*delta,'y<')
title('Quadratic regression for N. Atlantic ocean')
hold on
%figure()
%plot(TEMP_NATL(ido),NF_NATL(ido),'cy<')
%hold on
[p,ErrorEst] = polyfit(TEMP_NATL(ido),NF_NATL(ido),2);
[NF_fit,delta] = polyval(p,TEMP_NATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_NATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

% ======================================================================
fprintf('S. Atlantic ocean regression \n')
% select data only in the S. Atlantic ocean for regression. 
ido = find(NF_SATL~=0);
RG_SATL = fitlm(TEMP_SATL(ido),NF_SATL(ido));
fprintf('R^2 = %0.3f; \n',  RG_SATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SATL.Coefficients.pValue(2))
[p,ErrorEst] = polyfit(TEMP_SATL(ido),NF_SATL(ido),1);
[NF_fit,delta] = polyval(p,TEMP_SATL(ido),ErrorEst);

figure()
plot(TEMP_SATL(ido),NF_SATL(ido),'k+',...
  TEMP_SATL(ido),NF_fit,'r*',...
  TEMP_SATL(ido),NF_fit+2*delta,'y<',...
  TEMP_SATL(ido),NF_fit-2*delta,'y<')
title('linear regression for S. Atlantic ocean')

%plot(TEMP_SATL(ido),NF_SATL(ido),'g*')
%hold off
[p,ErrorEst] = polyfit(TEMP_SATL(ido),NF_SATL(ido),2);
[NF_fit,delta] = polyval(p,TEMP_SATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_SATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

% =======================================================================
fprintf('N. Pacific ocean regression \n')
% select data only in the N. Pacific ocean for regression. 
ido = find(NF_NPAC~=0);
RG_NPAC = fitlm(TEMP_NPAC(ido),NF_NPAC(ido));
fprintf('R^2 = %0.3f; \n',  RG_NPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NPAC.Coefficients.pValue(2))
[p,ErrorEst] = polyfit(TEMP_NPAC(ido),NF_NPAC(ido),1);
[NF_fit,delta] = polyval(p,TEMP_NPAC(ido),ErrorEst);

figure()
plot(TEMP_NPAC(ido),NF_NPAC(ido),'k+',...
  TEMP_NPAC(ido),NF_fit,'r*',...
  TEMP_NPAC(ido),NF_fit+2*delta,'y<',...
  TEMP_NPAC(ido),NF_fit-2*delta,'y<')
title('Quadratic regression for N. Pacific ocean')
hold on
[p,ErrorEst] = polyfit(TEMP_NPAC(ido),NF_NPAC(ido),2);
[NF_fit,delta] = polyval(p,TEMP_NPAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_NPAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)
%figure()
%plot(TEMP_NPAC(ido),NF_NPAC(ido),'cy<')
%hold on
% =======================================================================
fprintf('S. Pacific ocean regression \n')
% select data only in the S. Pacific ocean for regression. 
ido = find(NF_SPAC~=0);
RG_SPAC = fitlm(TEMP_SPAC(ido),NF_SPAC(ido));
fprintf('R^2 = %0.3f; \n',  RG_SPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SPAC.Coefficients.pValue(2))
%plot(TEMP_SPAC(ido),NF_SPAC(ido),'g*')
%hold off
[p,ErrorEst] = polyfit(TEMP_SPAC(ido),NF_SPAC(ido),1);
[NF_fit,delta] = polyval(p,TEMP_SPAC(ido),ErrorEst);

figure()
plot(TEMP_SPAC(ido),NF_SPAC(ido),'k+',...
  TEMP_SPAC(ido),NF_fit,'r*',...
  TEMP_SPAC(ido),NF_fit+2*delta,'y<',...
  TEMP_SPAC(ido),NF_fit-2*delta,'y<')
title('linear regression for S. Pacific ocean')
hold off
[p,ErrorEst] = polyfit(TEMP_SPAC(ido),NF_SPAC(ido),2);
[NF_fit,delta] = polyval(p,TEMP_SPAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_SPAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

