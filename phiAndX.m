function [h,g] = phiAndX(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorList,baseParams,nSamps)

    %Nicholas Szczecinski
    %West Virginia University
    %13 October 2020

    lPower = linspace(lPowerMin,lPowerMax,nSamps);
    tPower = linspace(tPowerMin,tPowerMax,nSamps);

    %Generate vectors for t and l
    l = 10.^lPower;
    t = 10.^tPower;

    %Generate 2D grids for t and l
    [T,L] = meshgrid(t,l);
    
    %Read in the parameter values when L=1.
    m0 = baseParams(1);
    c0 = baseParams(2);
    k0 = baseParams(3);
    s = baseParams(4);
    g = 10;
    
    %Because the torsional analogs of m, c, and k are complicated, define
    %them as their own function, which we can use consistently everywhere.
    mOfL = @(L) 1/3*m0*L.^5;
    cOfL = @(L) c0*s^2*L.^3;
    kOfL = @(L) k0*s^2*L.^3 + m0*g/2*L.^4;
    
    %To get a basic sampling of m, c, and k over the sample space, just
    %plug in L and save the output.
    m = mOfL(L);
    c = cOfL(L);
    k = kOfL(L);

    %Compute natural frequency and damping ratio over the T and L grid.
    omega_n = sqrt(k./m);
    zeta = c./(2*omega_n.*m);  

    %r is the ratio between forcing frequency and natural frequency. This
    %is effectively the Strouhal number.
    omega = 2*pi./T;
    r = omega./omega_n;
    
    %phi is the phase angle between the forcing function and the appendage
    %displacement. 
    phi = atan2d(2*zeta.*r,1-r.^2);
    
    %X is the ratio between the input force divided by the "static
    %deflection" force (i.e. normalized by k) and the displacement of the
    %appendage.
    X = 1./sqrt( (1 - r.^2).^2 + (2*zeta.*r).^2);
    
    %mag is the magnitude of the force required to move the appendage given
    %its properties k, c, and m and the forcing frequency omega.
    mag = sqrt( (k - m.*omega.^2).^2 + (c.*omega).^2 );
    
    %Define a function that computes the damping ratio zeta as a function
    %of L
    zetaOfL = @(L) cOfL(L)./(2*sqrt(kOfL(L).*mOfL(L)));
    
    %Numerically compute the length scale at which zeta=1, signifying the
    %boundary between underdamping and overdamping (i.e. critically
    %damped).
    LcritDamped = fzero(@(L)1-zetaOfL(L),[1e-3,1]);
    
    %Define a fine, logarithmically-spaced vector along the L dimension,
    %which we can use to draw boundaries between regions in a smooth, clean
    %manner.
    LforBnds = logspace(lPowerMin,lPowerMax,1000);
    mLo = mOfL(LforBnds);
    cLo = cOfL(LforBnds); %#ok<NASGU>
    kLo = kOfL(LforBnds);

    %Compute the natural period boundary.
    Tn = 2*pi*sqrt(mLo./kLo);

    %Compute the boundary between quasi-static and viscous regimes. This
    %occurs when the phase angle is 60 degrees, which would divide the
    %range of possible phase lags, [0 180], into three equal regions.
    %However, to facilitate analysis, these boundaries compute where 
    %tan(phi) = +/- 2. This corresponds to a phi value of about 57 degrees
    %instead of exactly 60 (and 117 instead of 120).
    %The "low boundary" is between quasi-static and viscous, and the "high
    %boundary" is between viscous and kinetic.
    TloBnd =  (2*pi*cOfL(LforBnds) + 2*pi*sqrt(cOfL(LforBnds).^2 + 16*kOfL(LforBnds).*mOfL(LforBnds)))./(4*kOfL(LforBnds));
    ThiBnd = (-2*pi*cOfL(LforBnds) + 2*pi*sqrt(cOfL(LforBnds).^2 + 16*kOfL(LforBnds).*mOfL(LforBnds)))./(4*kOfL(LforBnds));

    %Plot a preliminary plot wherein the background is shaded by the
    %forcing regime (kinetic, viscous, or quasi-static); boundaries between
    %these regions are plotted; and animal behavior "bars" are plotted over
    %top.
    h = figure;
    surf(T,L,phi,'edgealpha',0);
    view(0,90)
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('\phi')
    c = colorbar;
    c.Limits = [0,180];
    c.Ticks = 0:90:180;
    c.TickLabels = {'0\circ (Q.S.)','90\circ (Vis.)','180\circ (Kin.)'};
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    hold on
    LnatFreq = LforBnds(LforBnds >= LcritDamped);
    TnatFreq = Tn(LforBnds >= LcritDamped);
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'k--','linewidth',1)
    text(min(TnatFreq),min(LnatFreq),180,'T = T_n (\phi = 90\circ)','VerticalAlignment','bottom','Rotation',45)
    Tlo = TloBnd(TloBnd>10^tPowerMin);
    Llo = LforBnds(TloBnd>10^tPowerMin);
    plot3(Tlo,Llo,180+zeros(size(Tlo)),'k','linewidth',1)
    text(min(Tlo),min(Llo),180,'T = T_{ev} (\phi \approx 63\circ)','VerticalAlignment','bottom','Rotation',90)
    Thi = ThiBnd(ThiBnd>10^tPowerMin);
    Lhi = LforBnds(ThiBnd>10^tPowerMin);
    plot3(Thi,Lhi,180+zeros(size(Thi)),'k','linewidth',1)
    text(min(Thi),min(Lhi),180,'T = T_{vi} (\phi \approx 117\circ)','VerticalAlignment','bottom','Rotation',atand(1/2))
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    
    M = contourc(t,l,zeta,[1,1]);
    plot3(M(1,2:end),M(2,2:end),180+zeros(size(t)),'w--');
    
    for i=1:size(behaviorList,1)
        plot3(behaviorList{i,2},[1,1]*behaviorList{i,3},[200,200],'r','linewidth',1)
        text(behaviorList{i,2}(2)*behaviorList{i,4},behaviorList{i,3}*behaviorList{i,5},200,behaviorList{i},'fontweight','bold','fontsize',10,'Color','red')
    end
    set(h,'Position',[480.2000e+000 153.8000e+000 586.4000e+000 479.2000e+000])

    text(10^1.8,1.5*M(2,end),200,'underdamped','fontweight','bold','fontsize',12,'color','black','HorizontalAlignment','right')
    text(10^1.8,1/1.5*M(2,end),200,'overdamped','fontweight','bold','fontsize',12,'color','black','HorizontalAlignment','right')
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Compute some approximate boundaries between our regions. These
    %approximations have the advantage of being simple formulations and
    %straight lines in the log-log space.
    %T - time axis values
    %L - length axis values
    %e - elastic (potential)
    %v - viscous
    %i - inertial (kinetic)
    Tev = pi*c0/k0*[1,1];
    Lev = [10^lPowerMin,LcritDamped];
    
    Lvi = LforBnds(LforBnds <= LcritDamped);
    Tvi = 4*pi*(1/3*m0)/(c0*s^2)*Lvi.^2;

    Tlo = TloBnd(TloBnd>10^tPowerMin);
    Llo = LforBnds(TloBnd>10^tPowerMin);

    Thi = ThiBnd(ThiBnd>10^tPowerMin);
    Lhi = LforBnds(ThiBnd>10^tPowerMin);

    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    
    figure
    lw = 1;
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('\phi')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    
    hold on
    LnatFreq = LforBnds(LforBnds >= LcritDamped);
    TnatFreq = Tn(LforBnds >= LcritDamped);
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'k--','linewidth',1)
    text(min(TnatFreq),min(LnatFreq),180,'T = T_n (\phi = 90\circ)','VerticalAlignment','bottom','Rotation',45)
    
    plot3(2*pi*sqrt(1/3*m0*l.^5./(k0*s^2*l.^3 + 1/2*m0*10*l.^4)), l, 180+zeros(size(l)),':')
    
    Tlo = TloBnd(TloBnd>10^tPowerMin);
    Llo = LforBnds(TloBnd>10^tPowerMin);
    plot3(Tlo,Llo,180+zeros(size(Tlo)),'k','linewidth',1)
    text(min(Tlo),min(Llo),180,'T = T_{ev} (\phi \approx 63\circ)','VerticalAlignment','bottom','Rotation',90)
    Thi = ThiBnd(ThiBnd>10^tPowerMin);
    Lhi = LforBnds(ThiBnd>10^tPowerMin);
    plot3(Thi,Lhi,180+zeros(size(Thi)),'k','linewidth',1)
    text(min(Thi),min(Lhi),180,'T = T_{vi} (\phi \approx 117\circ)','VerticalAlignment','bottom','Rotation',atand(1/2))
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'k','linewidth',lw)
    plot3(Tev,Lev,180+zeros(size(Tev)),'k--','linewidth',lw)
    plot3(Tvi,Lvi,180+zeros(size(Tvi)),'k-.','linewidth',lw)
    view(0,90)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Plot the figure of animal behaviors on top of the dynamic regimes.
    
    %Z is a remapping of zeta for color-coding. When zeta is 0, Z is 1.
    %When zeta approaches infinity, Z is -1. Z can be thought of as an
    %"underdamping metric". Z>0->underdamped, Z<0->overdamped.
    Z = 1-zeta./(1 + zeta);
    g = figure;
    
    surf(T,L,phi,'edgealpha',0,'facealpha',1)
    colormap autumn
    
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    hold on
    lw = 1;
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'k','linewidth',lw)
    plot3(Tev,Lev,180+zeros(size(Tev)),'k--','linewidth',lw)
    plot3(Tvi,Lvi,180+zeros(size(Tvi)),'k-.','linewidth',lw)
    plot3([10^tPowerMin,10^tPowerMax],LcritDamped*[1,1],[180,180],'k:','linewidth',lw)
    
    %Compute the limits due to bone and muscle stress.
    levs = [8,8]; %100 MPa because of log scaling
    boneStressContour = contour(T,L,log10(mag./(.012*L.^2)./(s*L)),levs);
    levs = [7,7]; %10 MPa because of log scaling  
    tendonStressContour = contour(T,L,log10(mag./(.012*L.^2)./(s*L)),levs);
    
    plot3(boneStressContour(1,2:end),boneStressContour(2,2:end),180+zeros(size(boneStressContour(1,2:end))),'--','color',.75*[1,1,1],'linewidth',2)
    plot3(tendonStressContour(1,2:end),tendonStressContour(2,2:end),180+zeros(size(tendonStressContour(1,2:end))),':','color',.75*[1,1,1],'linewidth',2)
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic                viscous                quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    
    for i=1:size(behaviorList,1)
        plot3(behaviorList{i,2},[1,1]*behaviorList{i,3},[200,200],'k','linewidth',1.5)
        text(behaviorList{i,2}(2)*behaviorList{i,4},behaviorList{i,3}*behaviorList{i,5},200,behaviorList{i},'fontweight','bold','fontsize',8,'Color','black','fontname','Arial')
    end
    set(g,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])

    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-2.3,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    drawnow


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Define time and length scales of interest, which we will label in our
    %T and L plot with the letters A, B, C, D, and E.
    Lletter(1) = 1e-2;
    Tletter(1) = 3;
    
    Lletter(2) = 1.7;
    Tletter(2) = 3;
    
    Lletter(3) = 1.7;
    Tletter(3) = 1.5; %0.5;

    Lletter(4) = 2e-3;
    Tletter(4) = 1e-2;

    Lletter(5) = .03;
    Tletter(5) = 10^-2.5;
    
    letterCell = {'D','C','B','E','F'};
    letterColor = {'black','black','black','black','black'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    g3 = figure;
    spSurf = subplot(3,3,[1,2,4,5]);
    surf(T,L,phi,'edgealpha',0,'facealpha',1)
    colormap autumn
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('A')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    % surfAx.ZScale = 'log';
    hold on
    lw = 1;
    bndColor = 0*[1,1,1];
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'color',bndColor,'linewidth',lw,'linestyle','-')
    plot3(Tev,Lev,180+zeros(size(Tev)),'color',bndColor,'linewidth',lw,'linestyle','--')
    plot3(Tvi,Lvi,180+zeros(size(Tvi)),'color',bndColor,'linewidth',lw,'linestyle','-.')
    plot3([10^tPowerMin,10^tPowerMax],LcritDamped*[1,1],[180,180],'color',bndColor,'linewidth',lw,'linestyle',':')
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic          viscous           quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    cbar.Position(1) = .53;

    set(g3,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])
    
    surfAx.XTick = logspace(tPower(1),tPower(end),tPower(end)-tPower(1)+1);
    surfAx.YTick = logspace(lPower(1),lPower(end),lPower(end)-lPower(1)+1);
    surfAx.Position(4) = surfAx.Position(4) - 0.05;
    surfAx.Position(2) = surfAx.Position(2) + 0.05;
    cbar.Position(2) = surfAx.Position(2);
    cbar.Position(4) = surfAx.Position(4);
    
    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-1.8,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    
    sp1 = [];
    sp2 = [];
%     %QS, overdamped
%     [tSim{1},xSim{1},~,~,~,~,~,Fapp{1},~] = simulateJointResponse(k0,c0,m0,s,Lletter(1),3,Tletter(1),15*Tletter(1),sp1,sp2);
%     title(letterCell{1})
%     %QS, underdamped
%     [tSim{2},xSim{2},~,~,~,~,~,Fapp{2},~] = simulateJointResponse(k0,c0,m0,s,Lletter(2),3,Tletter(2),15*Tletter(2),sp1,sp2);
%     title(letterCell{2})
%     %Inertial, underdamped
%     [tSim{3},xSim{3},~,~,~,~,~,Fapp{3},~] = simulateJointResponse(k0,c0,m0,s,Lletter(3),3,Tletter(3),15*Tletter(3),sp1,sp2);
%     title(letterCell{3})
%     %viscous, overdamped
%     [tSim{4},xSim{4},~,~,~,~,~,Fapp{4},~] = simulateJointResponse(k0,c0,m0,s,Lletter(4),3,Tletter(4),15*Tletter(4),sp1,sp2);
%     title(letterCell{4})
%     %Inertial, overdamped
%     [tSim{5},xSim{5},~,~,~,~,~,Fapp{5},~] = simulateJointResponse(k0,c0,m0,s,Lletter(5),3,Tletter(5),15*Tletter(5),sp1,sp2);
%     title(letterCell{5})
    
    for i=1:5
        [tSim{i},xSim{i},~,~,~,~,~,Fapp{i},~] = simulateJointResponse(k0,c0,m0,s,Lletter(i),3,Tletter(i),15*Tletter(i),sp1,sp2);
        title(letterCell{i})
    end
    
    spInd = [9,6,3,8,7];
    
    figure(g3)
    
    for i=1:5
        %Label on the T vs. L plot
        subplot(spSurf)
        hold on
        box on
        text(Tletter(i),Lletter(i),200,letterCell{i},'color',letterColor{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold')
        
        %Pick the subplot
        spWL = subplot(3,3,spInd(i));
        
        %Plot work loop
        title(letterCell{i})
        pts = (tSim{i} <= Tletter(i) & Fapp{i}(tSim{i}) > 0);
        hold on
        
        cvec = autumn(3);
        %x sim is negative because of how a muscle would shorten. Consider:
        %extensor causes positive torque. But when the joint is moving in
        %the positive direction, the muscle is shortening (negative
        %direction). Thus the joint angle should have the opposite sign to
        %measure the work done by the actuator.   
        
        t = tSim{i}(pts);
        F = max(Fapp{i}(t),0);
        x = -xSim{i}(pts);
        
        [maxF,ind] = max(F);
        if x(ind) < -0.1
            %quasi-static
        
            %Find the index at which the minimum position is reached. Split
            %up F and x before and after this point.
            [~,ind] = min(x);
            ftop = F(1:ind);
            xtop = x(1:ind);
            fbot = F(ind:end);
            xbot = x(ind:end);

            %Sample the "top" and "bottom" of the loop at the same grid,
            %allowing us to use the area() plot.
            Xbot = -1:.01:1;
            Fbot = interp1(xbot,fbot,Xbot);
            Fbot(isnan(Fbot)) = 0;
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Ftop)) = 0;

            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == -1) = max(F);
            Ftop(Xtop == -1) = max(F);
            
            %Plot the area under these curves and color code them based on
            %where the energy goes.
            area(Xbot',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            area(Xbot',Fbot','edgealpha',0,'facecolor',cvec(1,:));
            plot(Xbot,Fbot,'k','linewidth',1)
            plot(Xtop,Ftop,'k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/2);
            botStrt = [xbot(botArrLoc),fbot(botArrLoc)];
            botStop = [xbot(botArrLoc+1),fbot(botArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
            
            topArrLoc = ceil(length(xtop)/3);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        elseif x(ind) > 0.1
            %inertial (kinetic)
            
            %Find the index at which the minimum position is reached. Split
            %up F and x before and after this point.
            [~,ind] = max(x);
            ftop = F(ind:end);
            xtop = x(ind:end);
            fbot = F(1:ind);
            xbot = x(1:ind);

            %Sample the "top" and "bottom" of the loop at the same grid,
            %allowing us to use the area() plot.
            Xbot = -1:.01:1;
            Fbot = interp1(xbot,fbot,Xbot);
            Fbot(isnan(Fbot)) = 0;
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Ftop)) = 0;
            
            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == 1) = max(F);
            Ftop(Xtop == 1) = max(F);

            %Plot the area under these curves and color code them based on
            %where the energy goes.
            area(Xbot',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            area(Xbot',Fbot','edgealpha',0,'facecolor',cvec(3,:));
            plot(Xbot',Ftop','k','linewidth',1)
            plot(Xbot',Fbot','k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/3);
            botStrt = [xbot(botArrLoc),fbot(botArrLoc)];
            botStop = [xbot(botArrLoc+1),fbot(botArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
            
            topArrLoc = ceil(length(xtop)/2);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        else
            %viscous
            
            [~,ind] = min(x);
            ftop = F(1:ind);
            xtop = x(1:ind);
            fbot = F(ind:end);
            xbot = x(ind:end);

            Xbot = -1:.01:1;
            Fbot = zeros(size(Xbot));
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Xtop)) = 0;
            
            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == 1) = 0;
            Ftop(Xtop == 1) = 0;
            Fbot(Xbot == -1) = 0;
            Ftop(Xbot == -1) = 0;

            area(Xtop',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            plot(Xtop',Ftop','k','linewidth',1)
            plot(Xbot',Fbot','k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/2);
            botStrt = [-1,0];
            botStop = [.1,0];
            arrObj = arrow(botStrt,botStop);
            arrow fixlimits
            
            topArrLoc = ceil(length(xtop)/2);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        end        

        if ~strcmp(letterCell{i},'B') && ~strcmp(letterCell{i},'C')
            xlabel('Angle, \theta (rad)')
        end
        if ~strcmp(letterCell{i},'D') && ~strcmp(letterCell{i},'E')
            ylabel('Moment, M (Nm)')
        end
        
        title(letterCell{i})
        
        box on
        set(gca, 'Layer', 'Top')
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    g4 = figure;
    sp = subplot(3,3,[1,2,4,5]);
    surf(T,L,phi,'edgealpha',0,'facealpha',1)
    colormap autumn
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('A')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    hold on
    lw = 1;
    bndColor = 0*[1,1,1];
    plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'color',bndColor,'linewidth',lw,'linestyle','-')
    plot3(Tev,Lev,180+zeros(size(Tev)),'color',bndColor,'linewidth',lw,'linestyle','--')
    plot3(Tvi,Lvi,180+zeros(size(Tvi)),'color',bndColor,'linewidth',lw,'linestyle','-.')
%     plot3(TloBnd,LforBnds,180+zeros(size(LforBnds)),'color',bndColor,'linewidth',lw)
%     plot3(ThiBnd,LforBnds,180+zeros(size(LforBnds)),'color',bndColor,'linewidth',lw)
    
%     keyboard
    

    plot3([10^tPowerMin,10^tPowerMax],LcritDamped*[1,1],[180,180],'color',bndColor,'linewidth',lw,'linestyle',':')
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic              viscous       quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    cbar.Position(1) = .53;
    

    set(g4,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])
    surfAx.Position(4) = surfAx.Position(4) - 0.05;
    surfAx.Position(2) = surfAx.Position(2) + 0.05;
    cbar.Position(2) = surfAx.Position(2);
    cbar.Position(4) = surfAx.Position(4);
    
    surfAx.XTick = logspace(tPower(1),tPower(end),tPower(end)-tPower(1)+1);
    surfAx.YTick = logspace(lPower(1),lPower(end),lPower(end)-lPower(1)+1);

    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-1.8,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    stimMultiplier = [1;1;3;1;1]+.5;
    numCycles = [5;10;10;10;20];
    
    sp2 = [];
    
    for i=1:5
        sp1 = subplot(3,3,spInd(i));
        [tSim{i},xSim{i},~,~,~,~,~,Fapp{i},~] = simulateJointResponse(k0,c0,m0,s,Lletter(i),numCycles(i),Tletter(i),stimMultiplier(i)*Tletter(i),sp1,sp2);
        title(letterCell{i})
        
        if ~strcmp(letterCell{i},'B') && ~strcmp(letterCell{i},'C')
            xlabel('Time, t (n.d.)')
        else
            xlabel('')
        end
        if ~strcmp(letterCell{i},'D') && ~strcmp(letterCell{i},'E')
            ylabel('Angle, \theta (rad)')
        else
            ylabel('')
        end
    end
    
    %Add titles to subplot
    subplot(sp)
    for i=1:5
        text(Tletter(i),Lletter(i),200,letterCell{i},'color',letterColor{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold')
    end 
    
    
end