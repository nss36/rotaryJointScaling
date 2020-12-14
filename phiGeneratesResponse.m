function [h,g] = phiGeneratesResponse(phi,x,t)
    
    %Nicholas Szczecinski
    %West Virginia University
    %13 October 2020

    %Start a 3D plot of the forces over time. Plot the response of each
    %system (based on its value for phi) in a staggered way.
    h = figure;
    plot3(zeros(size(t)),t,x)
    hold on
    grid on
    xlabel('\phi')
    ylabel('t')
    zlabel('F')
    ax = gca;
    ax.XTick = phi;
    view(50,30);
    xlim([min(phi),max(phi)])
    ylim([min(t),max(t)])
    
    %Start a 2D fgure with subplots of the forces over time. Plot the
    %response of each system on its own axes.
    g = figure;
    
    xdot = centeredDiff(t,x);
    xddot = centeredDiff(t,xdot);

    %To help the user debug, we will plot the calculated values of k, c,
    %and m. Plot table column headings.
    fprintf('phi\tk\t\tc\t\tm\n');
    
    %We wish to color code each trace by phi. In order to perform a
    %discrete interpolation of the colors, 1) define a colormap with 181
    %values (phi = 0 to phi = 180 in 1 degree increments); 2) for each phi
    %value, find the nearest integer degree; 3) use this value to map the
    %color.
    allLineColors = autumn(181);
    lineColors = NaN(length(phi),3);
    for i=1:length(phi)
        lineColors(i,:) = allLineColors(round(phi(i))+1,:);
    end
    
    Fmax = 0;
    
    %For a given phi value, there are infinite combinations of m, c, and k.
    %However, we can arbitrarily pick one such combination based on the
    %vector definition of phi.
    for i=1:length(phi)
        %Find the real and imaginary components of the phase angle between
        %the forcing function and the limb displacement. This angle depends
        %on the ratios between k, c, and m.
        re = cosd(phi(i));
        im = sind(phi(i));
        
        %Since re exists [-1,1], we can ensure m is nonnegative by setting it
        %equal to 1-re. Also, k-m=re, according to the x-axis vector sum of
        %forces, so we define k and m in the following way:
        k = 1;
        m = 1-re;
        c = im;
        
        %Add these values to the table printed out to the display.
        fprintf('%i\t%2.2f\t%2.2f\t%2.2f\n',phi(i),k,c,m);
              
        %Calculate the total force as a function of the passive forces.
        F = k*x + c*xdot + m*xddot;
                
        %In order to scale axes later, keep track of the maximum force. If
        %the maximum force in this loop is higher than it's ever been,
        %overwrite the previous max force.
        if max(F) > Fmax
            Fmax = max(F);
        end

        %Plot these forces on the 3D plot. Plot the agonist force as a
        %solid trace, plot the antagonist (rectified negative) force as a
        %dotted trace.
        figure(h)
        hold on
        pp1 = plot3(phi(i)+zeros(size(t)),t,max(F,0),'linewidth',2,'color',lineColors(i,:)); %#ok<NASGU>
        if any(F < 0)
            pp2 = plot3(phi(i)+zeros(size(t)),t,max(-F,0),':','linewidth',2,'color',lineColors(i,:)); %#ok<NASGU>
        end
        
        %Plot these forces on the 2D subplot.
        figure(g)
        hold on
        subplot(length(phi)+1,1,i+1)
        pp1 = plot(t,max(F,0),'linewidth',2,'color',lineColors(i,:)); %#ok<NASGU>
        if any(F < 0)
            hold on
            pp2 = plot(t,max(-F,0),':','linewidth',2,'color',lineColors(i,:)); %#ok<NASGU>
        end
        
    end
    
    %Adjust the axes, scaling, and labels on the subplots.
    figure(g)
    hold on
    subplot(length(phi)+1,1,1)
    plot(t,x,'k')
    ylabel('x')
    hold on
    for i=1:length(phi)
        subplot(length(phi)+1,1,i+1)
        ylim([0,Fmax]);
        ylabel('F')
        title(['\phi = ',num2str(phi(i)),'^o'])
        
        if i == length(phi)
            xlabel('t')
        end
    end
   
end