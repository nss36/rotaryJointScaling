function ff = NEURONEXscalingPlot(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorList,baseParams,nSamps)

    %Nicholas Szczecinski
    %West Virginia University
    %13 October 2020

    %Vectors for sampling the time- and length-scales
    t = logspace(tPowerMin,tPowerMax,nSamps);
    l = logspace(lPowerMin,lPowerMax,nSamps);

    
    %Grid up T and L.
    [T,L] = meshgrid(t,l);
    
    m0 = baseParams(1);
    c0 = baseParams(2);
    k0 = baseParams(3);
    s = baseParams(4);
    g = 10;

    %If m scales L^3 and k scales L, then the natural frequency scales 1/L.
    %Therefore, the natural period scales L. These data will be used to plot
    %the boundary between inertial and quasi-static regimes.
    mOfL = 1/3*m0*l.^5;
    kOfL = k0*s^2*l.^3 + m0*g/2*l.^4;
    cOfL = c0*s^2*l.^3;

    tNatural = 2*pi*sqrt(mOfL./kOfL);
    lNatural = l;
    
    %Compute the damped natural frequency. 
    zetaOfL = cOfL./(2*sqrt(kOfL).*sqrt(mOfL));
    tDampedNatural = tNatural./sqrt(1 - zetaOfL.^2);
    
    %Remove points that will be out of bounds.
    lNatural(tNatural < min(t)) = [];
    tNatural(tNatural < min(t)) = [];
    
    lDampedNatural = l;
    tDampedNatural(abs(imag(tDampedNatural)) > 0) = NaN;

    %Remove NaNs from vectors representing the damped natural period.
    tDampedNatural(isnan(lDampedNatural)) = [];
    lDampedNatural(isnan(lDampedNatural)) = [];
    lDampedNatural(isnan(tDampedNatural)) = [];
    tDampedNatural(isnan(tDampedNatural)) = [];

    %Do some other filtering that makes plotting easier.
    [lMin,minInd] = min(lDampedNatural);
    tMin = tDampedNatural(minInd);

    if minInd == 1
        tDampedNatural = [max(t),tDampedNatural];
        lDampedNatural = [lMin,lDampedNatural];
    else
        tDampedNatural = [tDampedNatural,max(t)];
        lDampedNatural = [lDampedNatural,lMin];
    end
    
    %%% START PLOTTING %%%
    %Heat map. This was a crowd pleaser while assembling the
    %NeuroNex.

    %Compute how the effective m, c, and k scale with L, and how each force
    %scales with L and T. These will be necessary for plotting the heatmap.
    M = 1/3*m0.*L.^5;
    C = c0*s^2.*L.^3;
    K = k0*s^2.*L.^3 + m0*g/2.*L.^4;
    Z = C./(2*sqrt(K).*sqrt(M));
    MA = M./T.^2; %inertial torque scaling
    KX = K; %restoring torque scaling
    
    %Ratio of inertial over elastic force (iOE)
    iOE = real(log10(MA./KX) .* (abs(imag(MA./KX)) <= 0));
    iOEnorm = (iOE - min(iOE(:)))./range(iOE(:));

    %Amount by which the damped natural frequency deviates from the
    %undamped natural frequency (if overdamped, this is infinite).
    dampedNaturalDeviation = M./K.*(1 - 1./sqrt(1 - Z.^2));
    
    viscousDomain = real(log10(dampedNaturalDeviation) .* (abs(imag(dampedNaturalDeviation)) <= 0));
    vDnorm = (viscousDomain - min(viscousDomain(:)))/range(viscousDomain(:));
    
    if all(isnan(vDnorm))
        %The entire domain being tested is overdamped.
      
        keyboard
    end
    
    %Using the y cb cr color scheme worked nicely because it is a 2D color
    %scheme, and mapped these values naturally to colors: https://en.wikipedia.org/wiki/YCbCr
    color = ycbcr2rgb(iOEnorm,vDnorm);

    ff = figure;
    subplot(3,3,[1,2,4,5,7,8])
    surf(T,L,zeros(size(T)),color,'edgealpha',0)
    hold on
    view(0,90)
    %Line styles are totally up for debate. There was concern that too harsh of
    %lines would suggest rigid boundaries that don't exist. On the other hand,
    %the point is to categorize! So feel free to experiment.
    plot(tNatural,lNatural,'k-.','linewidth',1)
    plot(tDampedNatural,lDampedNatural,'k:','linewidth',1)
    xlabel('cycle time (s)','fontsize',8)
    ylabel('length scale (L \propto 1/\omega_n)','fontsize',8)
    fa = gca;
    fa.XScale = 'log';
    fa.YScale = 'log';
    fa.XAxis.FontSize = 8;
    fa.YAxis.FontSize = 8;
    
    title('Body size and behavioral speed determine dynamic scale','fontsize',10)
    for i=1:size(behaviorList,1)
    
        %v3: new format, plus drawing lines instead of single markers.
        plot(behaviorList{i,2},[1,1]*behaviorList{i,3},'k','linewidth',1,'markersize',8)
        text(behaviorList{i,2}(2)*behaviorList{i,4},behaviorList{i,3}*behaviorList{i,5},behaviorList{i,1},'fontsize',8)
    
    end
    
    text(5e-3,3e0,'Inertial','color',[1,1,1])
    text(1e0,1.5e-1,'Quasi-static','color',[1,1,1])
    text(5e-3,5e-4,'Viscous','color',[1,1,1])
    xlim([min(t),max(t)])
    ylim([min(l),max(l)])
    box on

    subplot(3,3,[6,9])
    hold on
    [X,Y] = meshgrid(0:.01:1,0:.01:1);
    colorKey = ycbcr2rgb(X,Y);
    surf(X,Y,zeros(size(X)),colorKey,'edgealpha',0)
    title('Key','fontsize',8)
    fk = gca;
    fk.XDir = 'reverse';
    fk.YDir = 'reverse';
    xLabelList = [.75,.25,.5];
    yLabelList = [.2,.2,.8];
    labelList = {{'inertia','dominates'},{'elasticity','dominates'},{'viscosity','dominates'}};
    for i=1:3
        text(xLabelList(i),yLabelList(i),labelList{i},'fontsize',8,'HorizontalAlignment','center','color',[1,1,1])
    end
    fk.Position(2) = .2;
    fk.Position(4) = .4;
    fk.XAxis.Visible = 'off';
    fk.YAxis.Visible = 'off';
    xticklabels({})
    yticklabels({})
    ytickangle(90)
    view(0,90)
    box off

    set(ff,'Position',[551 200 622 363])
    set(ff,'renderer','Painters')
    
    
end

