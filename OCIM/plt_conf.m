function plt_conf(x,y,txt,islog,SB)

% this function is used to do second order regression 
% and plot model with observation curve with 95% conficence 
% interval. x is environment parameter, and y is observation.
% if the parameter that has been taken log, let islog = 1.
% SB is a structure that contains Basion specific data.
if nargin >4
        x = [SB.NATL_x',SB.SATL_x',SB.NPAC_x',SB.SPAC_x',SB.IND_x'];
        y = [SB.NATL_y',SB.SATL_y',SB.NPAC_y',SB.SPAC_y',SB.IND_y'];
    end

    [p,ErrorEst] = polyfit(x,y,2);
    [y_fit,delta] = polyval(p,x,ErrorEst);
    [r2,rms] = rsquare(y,y_fit);

    % Compute the real roots and determine the extent of the data.
    r      = roots(p)'; 		% Roots as a row vector.
    real_r = r(imag(r) == 0); 	% Real roots.
    alpha  = 0.05;              % Significant level

    % assure that the data are row vectors.
    x = reshape(x,1,length(x));
    y = reshape(y,1,length(y));

    % Extent of the data.
    mx = min(x);
    Mx = max(x);
    my = min(y);
    My = max(y);

    % Scale factors for plotting.
    sx = 0.05*(Mx-mx);
    sy = 0.05*(My-my);

    % Plot the data, the fit
    if nargin>4
        hdata1 = plot(SB.NATL_x,SB.NATL_y,'gd','markersize',3);
        hold on
        hdata2 = plot(SB.SATL_x,SB.SATL_y,'kd','markersize',3);
        hold on
        hdata3 = plot(SB.NPAC_x,SB.NPAC_y,'c*','markersize',3);
        hold on
        hdata4 = plot(SB.SPAC_x,SB.SPAC_y,'r*','markersize',3);
        hold on
        hdata5 = plot(SB.IND_x,SB.IND_y,'b+','markersize',3);
    else 
        hdata = plot(x,y,'md','MarkerSize',3,'LineWidth',2);
    end 
    hold on
    xfit = mx-sx:0.01:Mx+sx;
    yfit = polyval(p,xfit);
    hfit = plot(xfit,yfit,'b-','LineWidth',2);

    % Add prediction intervals to the plot.
    [Y,DELTA] = polyconf(p,xfit,ErrorEst,'alpha',alpha);
    hconf = plot(xfit,Y+DELTA,'b--');
    plot(xfit,Y-DELTA,'b--')

    % Display the polynomial fit and the real roots.
    approx_p = round(100*p)/100; % Round for display.

    % add regression funtion
    txt_title = texlabel(sprintf('%s = %+2.2f%s^2%+2.2f%s%+2.2f',txt.target,...
    approx_p(1),txt.parm,approx_p(2),txt.parm,approx_p(3)));

    R2 = texlabel(sprintf('R^2 = %0.2f',r2));
    titlename = [txt.Basion ', ' txt_title, ' , ' R2];
    title(titlename,'fontsize',10)
    if islog == 1
        txt_xlabel = texlabel(sprintf('LOG_10 (%s), %s',txt.parm_longname,txt.xunit));
        txt_ylabel = texlabel(sprintf('LOG_10 (%s), %s',txt.target_longname,txt.yunit));
    else
        txt_xlabel = texlabel(sprintf('%s, %s',txt.parm_longname,txt.xunit));
        txt_ylabel = texlabel(sprintf('LOG_10 (%s), %s',txt.target_longname,txt.yunit));
    end
    xlabel(txt_xlabel);
    ylabel(txt_ylabel);
    % Add a legend.
    if nargin > 4
        leg1 = legend([hdata1,hdata2,hdata3,hdata4,hdata5],...
        'N. Atlantic','S. Atlantic','N. Pacific','S. Pacific','Indian');
        set(leg1,'fontsize',14,'box','off');
        ch1 = get(leg1, 'Children');
        set(ch1(:), 'MarkerSize', 16);      % Color only
    else 
        leg2 = legend([hdata,hfit,hconf],...
        'Data','Fit','95% Prediction Intervals');
        set(leg2,'fontsize',14,'box','off');
        ch2 = get(leg2, 'Children');
        set(ch2(:), 'MarkerSize', 16);      % Color only

    end 
