function h = scalingInverseDynamicsAcrossScales(tPowerMin,tPowerMax,lPowerMin,lPowerMax,baseParams,nSamps)

    %Plot the relative contribution of each force/torque for a given motion.
    %For each point on the T and L plot, plot the displacement and each
    %force over time. Scale the force to be on the same order as the
    %displacement (for visualization purposes).

    lPower = linspace(lPowerMin,lPowerMax,nSamps);
    tPower = linspace(tPowerMin,tPowerMax,nSamps);

    l = 10.^lPower;
    t = 10.^tPower;

    %Make gridded values for T and L.
    [T,L] = meshgrid(t,l);

    %Specify the shape of the limb displacement. This can be sinusoidal;
    %sinrise (1 - cos(omega*t)); expwave (wave constructed of exponential
    %rises and falls); or erf (error function rise; corresponds to Gaussian
    %velocity profile).
    motionType = 'sinusoidal'; 
    %motionType = 'sinrise';
    %motionType = 'erf';
    %motionType = 'exponential';
    %motionType = 'expwave';

    %Amplitude of motion
    D = 1;

    %Exponent of L that determines how each unit scales.
    kExp = 3;
    kGravExp = 4;
    cExp = 3;
    mExp = 5;
    
    m0 = baseParams(1);
    c0 = baseParams(2);
    k0 = baseParams(3);
    s = baseParams(4);
    g = 10;

    nrows = size(L,1);
    ncols = size(L,2);

    m = 1/3*m0*L.^mExp;
    c = c0*s^2*L.^cExp;
    k = k0*s^2*L.^kExp + m0*g/2*L.^kGravExp;

    %Assuming sinusoidal motion, calculate the nondimensional displacement
    %of the system (i.e. amplitude of motion, divided by (force amplitude
    %divided by spring stiffness)).
    omega = sqrt(k./m);
    zeta = c./(2*omega.*m);  

    r = 2*pi./omega./T;
    phi = atan2d(2*zeta.*r,1-r.^2);
    X = 1./sqrt( (1 - r.^2).^2 + (2*zeta.*r).^2);

    %Plot timecourses on a grid. This might not be great for manuscript
    %presentation, but it is useful as an exploratory tool.
    h = figure;
    set(h,'Position',[551 250 622 363])

    for i=1:nrows

        for j=1:ncols
            tVec = linspace(0,T(i,j),1000)';
            if strcmp(motionType,'exponential')
                %Exponential that starts at T(i,j)/10, and has a time
                %constant of 3*T(i,j)/10.
                xDes = D*(1 - exp(-max(tVec - T(i,j)/10,0)/(.3*T(i,j))));
            elseif strcmp(motionType,'sinusoidal')
                %Sin wave with the frequency 2*pi/T(i,j).
                omegaStim = 2*pi/T(i,j);
                xDes = D*sin(omegaStim*tVec);
            elseif strcmp(motionType,'sinrise')
                %One half-cycle of a sin wave with the frequency
                %2*pi/T(i,j).                
                xDes = D + zeros(size(tVec));
                omegaStim = 2*pi/T(i,j);
                xRise = D/2*(1-cos(omegaStim*tVec(tVec <= T(i,j)/2)));
                len = length(xRise);
                xDes(1:len/2) = 0;
                xDes((1+len/2):(1.5*len)) = xRise;
                xDes(1.5*len+1:end) = D;
            elseif strcmp(motionType,'erf')
                %Erf function motion (Gaussian velocity)
                omegaStim = 2*pi/T(i,j);
                xDes = D/2*(1+erf(sqrt(pi)/2*omegaStim*(tVec-T(i,j)/2)));
            elseif strcmp(motionType,'expwave')
                
                dt = tVec(2) - tVec(1);
                omegaStim = 3*(2*pi/T(i,j));
                u = D*2*(-1/2 + (double(sin(omegaStim*tVec) > 0)));
                tau = .1*2*pi/omegaStim;
                
                if dt > tau/3
                    error('dt is too coarse for the low pass filter.')
                end

                %Generate the exponential wave by smoothing the square wave
                %with a low-pass filter.
                xDes = zeros(size(tVec));
                xDes(1) = u(1);
                for n=2:length(u)
                    xDes(n) = xDes(n-1) + dt/tau*(u(n) - xDes(n-1));
                end
            else
                xDes = zeros(size(tVec));
            end

            %Compute the velocity and acceleration of the desired motion.
            xDotDes = centeredDiff(tVec,xDes);

            xDotDotDes = centeredDiff(tVec,xDotDes);

            %Inertial force, F = ma
            Finertial = m(i,j)*xDotDotDes;

            %Passive viscous forces, F = cv
            Fviscous = c(i,j)*xDotDes;

            %Passive elastic and gravitational forces, F = kx
            Felastic = k(i,j)*xDes;

            %Normalize the force to its maximum value for plotting.
            F = [Finertial,Fviscous,Felastic];
            F = F/max(F,[],'all');

            subplot(nrows,ncols,(nrows-i)*ncols + j); 
            hold on

            %Make a custom color list.
            lineColors = lines(5);
            lineColors(3,:) = [];
            tempPurp = lineColors(3,:);
            lineColors(3,:) = lineColors(4,:);
            lineColors(4,:) = tempPurp;

            %Plot each force over time.
            for n=1:3
                plot(tVec,F(:,n),'color',lineColors(n,:),'linewidth',1)
            end
            %Plot the sum of forces over time.
            plot(tVec,sum(F,2),'color',lineColors(4,:),'linestyle',':','linewidth',1)
            
            %Plot the displacement over time.
            plot(tVec,xDes,'k--','linewidth',1)
            grid on

            %Add a title with basic info about T, L, and the amplitude of
            %the force.
            title({sprintf('log(L) = %2.2f, log(T) = %2.2f',log10(L(i,j)),log10(T(i,j))),sprintf('log(X) = %2.2f, phi = %2.2f',log10(X(i,j)),phi(i,j))})
            
            xlim([0,max(tVec)])
            xlabel('t (s)')
            box on

            %Add just one figure legend, since it is the same everywhere.
            if i == nrows && j == ncols
                legend('inertial force','viscous force','elastic + grav. force','net force','displacement','location','southeast')
            end

        end

    end
end