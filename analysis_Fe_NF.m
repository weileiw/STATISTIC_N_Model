clc 
clear all
close all

load NF2d_OPT_N2P_tanh.mat
load DFe_from_Ben.mat
load ~/MOCM/DATA/transport_v4.mat
load ~/MOCM/DATA/M3d90x180x24v2.mat

DFe = DFe_3d(:,:,1);
ifix = find(NF2d(:)>0);

NF = log10(NF2d(ifix));
NF_z_score = (NF-mean(NF))/std(NF);

DFe = DFe(ifix);
DFe_z_score = (DFe-mean(DFe))/std(DFe);

ATL_srf = MSKS.ATL(:,:,1);
PAC_srf = MSKS.PAC(:,:,1);
IND_srf = MSKS.IND(:,:,1);
MED_srf = MSKS.MED(:,:,1);
ARC_srf = MSKS.ARC(:,:,1);

NF = NF_z_score;
DFe = DFe_z_score;

% apply masks to fixation.
NF_ATL = NF.*ATL_srf(ifix);
DFe_ATL = DFe.*ATL_srf(ifix);

NF_PAC = NF.*PAC_srf(ifix);
DFe_PAC = DFe.*PAC_srf(ifix);

NF_IND = NF.*IND_srf(ifix);
DFe_IND = DFe.*IND_srf(ifix);

NF_ARC = NF.*ARC_srf(ifix);
DFe_ARC = DFe.*ARC_srf(ifix);

% ==========================================================================
fprintf('global ocean regression \n')
% linear regression 
RG_GLOB = fitlm(DFe,NF,'linear');
fprintf('R^2 = %0.3f; \n',  RG_GLOB.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_GLOB.Coefficients.pValue(2))
% linear regression
[p,ErrorEst] = polyfit(DFe,NF,1);
[NF_fit,delta] = polyval(p,DFe,ErrorEst);

figure()
plot(DFe,NF,'+',...
  DFe,NF_fit,'g*',...
  DFe,NF_fit+2*delta,'r<',...
  DFe,NF_fit-2*delta,'r<')
title('linear regression for global ocean')

% quadratic regression
[p,ErrorEst] = polyfit(DFe,NF,2);
[NF_fit,delta] = polyval(p,DFe,ErrorEst);
[r2,rms] = rsquare(NF,NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(DFe,NF,'+',...
  DFe,NF_fit,'g*',...
  DFe,NF_fit+2*delta,'r<',...
  DFe,NF_fit-2*delta,'r<')
title('Quadratic regression for global ocean')

% ========================================================================
fprintf('Atlantic ocean regression \n')
% select data only in the Atlantic ocean for regression. 
ido = find(NF_ATL~=0);
RG_ATL = fitlm(DFe_ATL(ido),NF_ATL(ido));

fprintf('R^2 = %0.3f; \n',  RG_ATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_ATL.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(DFe_ATL(ido),NF_ATL(ido),2);
[NF_fit,delta] = polyval(p,DFe_ATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_ATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(DFe_ATL(ido),NF_ATL(ido),'+',...
  DFe_ATL(ido),NF_fit,'g*',...
  DFe_ATL(ido),NF_fit+2*delta,'r<',...
  DFe_ATL(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Atlantic ocean')

%==========================================================================
fprintf('Pacific ocean regression \n')
% select data only in the Pacific ocean for regression. 
ido = find(NF_PAC~=0);
RG_PAC = fitlm(DFe_PAC(ido),NF_PAC(ido));

fprintf('R^2 = %0.3f; \n',  RG_PAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_PAC.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(DFe_PAC(ido),NF_PAC(ido),2);
[NF_fit,delta] = polyval(p,DFe_PAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_PAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(DFe_PAC(ido),NF_PAC(ido),'+',...
  DFe_PAC(ido),NF_fit,'g*',...
  DFe_PAC(ido),NF_fit+2*delta,'r<',...
  DFe_PAC(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Pacific ocean')

%===========================================================================
fprintf('Indian ocean regression \n')
% select data only in the Indian ocean for regression. 
ido = find(NF_IND~=0);
RG_IND = fitlm(DFe_IND(ido),NF_IND(ido));

fprintf('R^2 = %0.3f; \n',  RG_IND.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_IND.Coefficients.pValue(2))

% quadratic regression
[p,ErrorEst] = polyfit(DFe_IND(ido),NF_IND(ido),2);
[NF_fit,delta] = polyval(p,DFe_IND(ido),ErrorEst);
[r2,rms] = rsquare(NF_IND(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

figure()
plot(DFe_IND(ido),NF_IND(ido),'+',...
  DFe_IND(ido),NF_fit,'g*',...
  DFe_IND(ido),NF_fit+2*delta,'r<',...
  DFe_IND(ido),NF_fit-2*delta,'r<')
title('Quadratic regression for the Indian ocean')

%===========================================================================
fprintf('Arctic ocean regression \n')
% select data only in the Arctic ocean for regression. 
ido = find(NF_ARC~=0);
RG_ARC = fitlm(DFe_ARC(ido),NF_ARC(ido));

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
SATL_srf = SATL(:,:,1);

NPAC_srf = NPAC(:,:,1);
SPAC_srf = SPAC(:,:,1);

% apply masks to fixation.
NF_NATL = NF.*NATL_srf(ifix);
DFe_NATL = DFe.*NATL_srf(ifix);
NF_SATL = NF.*SATL_srf(ifix);
DFe_SATL = DFe.*SATL_srf(ifix);

NF_NPAC = NF.*NPAC_srf(ifix);
DFe_NPAC = DFe.*NPAC_srf(ifix);
NF_SPAC = NF.*SPAC_srf(ifix);
DFe_SPAC = DFe.*SPAC_srf(ifix);

% =====================================================================
fprintf('N. Atlantic ocean regression \n')
% select data only in the N. Atlantic ocean for regression. 
ido = find(NF_NATL~=0);
RG_NATL = fitlm(DFe_NATL(ido),NF_NATL(ido));
fprintf('R^2 = %0.3f; \n',  RG_NATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NATL.Coefficients.pValue(2))

[p,ErrorEst] = polyfit(DFe_NATL(ido),NF_NATL(ido),1);
[NF_fit,delta] = polyval(p,DFe_NATL(ido),ErrorEst);

figure()
plot(DFe_NATL(ido),NF_NATL(ido),'k+',...
  DFe_NATL(ido),NF_fit,'r*',...
  DFe_NATL(ido),NF_fit+2*delta,'y<',...
  DFe_NATL(ido),NF_fit-2*delta,'y<')
title('Quadratic regression for N. Atlantic ocean')
hold on
%figure()
%plot(DFe_NATL(ido),NF_NATL(ido),'cy<')
%hold on
[p,ErrorEst] = polyfit(DFe_NATL(ido),NF_NATL(ido),2);
[NF_fit,delta] = polyval(p,DFe_NATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_NATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

% ======================================================================
fprintf('S. Atlantic ocean regression \n')
% select data only in the S. Atlantic ocean for regression. 
ido = find(NF_SATL~=0);
RG_SATL = fitlm(DFe_SATL(ido),NF_SATL(ido));
fprintf('R^2 = %0.3f; \n',  RG_SATL.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SATL.Coefficients.pValue(2))
[p,ErrorEst] = polyfit(DFe_SATL(ido),NF_SATL(ido),1);
[NF_fit,delta] = polyval(p,DFe_SATL(ido),ErrorEst);

figure()
plot(DFe_SATL(ido),NF_SATL(ido),'k+',...
  DFe_SATL(ido),NF_fit,'r*',...
  DFe_SATL(ido),NF_fit+2*delta,'y<',...
  DFe_SATL(ido),NF_fit-2*delta,'y<')
title('linear regression for S. Atlantic ocean')

%plot(DFe_SATL(ido),NF_SATL(ido),'g*')
%hold off
[p,ErrorEst] = polyfit(DFe_SATL(ido),NF_SATL(ido),2);
[NF_fit,delta] = polyval(p,DFe_SATL(ido),ErrorEst);
[r2,rms] = rsquare(NF_SATL(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

% =======================================================================
fprintf('N. Pacific ocean regression \n')
% select data only in the N. Pacific ocean for regression. 
ido = find(NF_NPAC~=0);
RG_NPAC = fitlm(DFe_NPAC(ido),NF_NPAC(ido));
fprintf('R^2 = %0.3f; \n',  RG_NPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_NPAC.Coefficients.pValue(2))
[p,ErrorEst] = polyfit(DFe_NPAC(ido),NF_NPAC(ido),1);
[NF_fit,delta] = polyval(p,DFe_NPAC(ido),ErrorEst);

figure()
plot(DFe_NPAC(ido),NF_NPAC(ido),'k+',...
  DFe_NPAC(ido),NF_fit,'r*',...
  DFe_NPAC(ido),NF_fit+2*delta,'y<',...
  DFe_NPAC(ido),NF_fit-2*delta,'y<')
title('Quadratic regression for N. Pacific ocean')
hold on
[p,ErrorEst] = polyfit(DFe_NPAC(ido),NF_NPAC(ido),2);
[NF_fit,delta] = polyval(p,DFe_NPAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_NPAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)
%figure()
%plot(DFe_NPAC(ido),NF_NPAC(ido),'cy<')
%hold on
% =======================================================================
fprintf('S. Pacific ocean regression \n')
% select data only in the S. Pacific ocean for regression. 
ido = find(NF_SPAC~=0);
RG_SPAC = fitlm(DFe_SPAC(ido),NF_SPAC(ido));
fprintf('R^2 = %0.3f; \n',  RG_SPAC.Rsquared.Adjusted)
fprintf('F statistic p-value = %5.5e \n',RG_SPAC.Coefficients.pValue(2))
%plot(DFe_SPAC(ido),NF_SPAC(ido),'g*')
%hold off
[p,ErrorEst] = polyfit(DFe_SPAC(ido),NF_SPAC(ido),1);
[NF_fit,delta] = polyval(p,DFe_SPAC(ido),ErrorEst);

figure()
plot(DFe_SPAC(ido),NF_SPAC(ido),'k+',...
  DFe_SPAC(ido),NF_fit,'r*',...
  DFe_SPAC(ido),NF_fit+2*delta,'y<',...
  DFe_SPAC(ido),NF_fit-2*delta,'y<')
title('linear regression for S. Pacific ocean')
hold off
[p,ErrorEst] = polyfit(DFe_SPAC(ido),NF_SPAC(ido),2);
[NF_fit,delta] = polyval(p,DFe_SPAC(ido),ErrorEst);
[r2,rms] = rsquare(NF_SPAC(ido),NF_fit);
fprintf('Quadratic regression R^2 = %0.3f; \n \n \n',  r2)

