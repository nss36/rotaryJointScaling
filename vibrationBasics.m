close all

r = logspace(-3,3,101);
z = logspace(-2,2,5);
t = 0:1e-3:30;
x0 = 1;

h = figure;

colors = prism(length(z)+1);
colors(1,:) = [];

spXofT = subplot(2,3,1:3);
xlabel('t')
ylabel('x')
title('A. Transient limb response timecourse')

spX = subplot(2,3,4);
ax = gca;
ax.XScale = 'log';
ylim([0,3])
xlabel('\omega/\omega_n')
ylabel('X')
xticks([1e-3,1e0,1e3])
title({'B. Steady state','response amplitude'})

spF = subplot(2,3,5);
ax = gca;
ax.XScale = 'log';
ylim([0,3])
xlim([1e-3,1e3])
xlabel('\omega/\omega_n')
ylabel('F_0')
xticks([1e-3,1e0,1e3])
title({'C. Steady state force required','for 1 rad amplitude'})

spPhi = subplot(2,3,6);
ax = gca;
ax.XScale = 'log';
ylim([0,180])
yticks(0:45:180)
xlabel('\omega/\omega_n')
ylabel('\phi (degrees)')
xticks([1e-3,1e0,1e3])
title({'D. Steady state','response phase angle'})

for i=1:length(z)
    
    subplot(spXofT)
    hold on
    
    if z(i) < 1
        A = 0;
        B = x0;
        plot(t,exp(-z(i)*t).*(A*sin(sqrt(1-z(i)^2)*t) + B*cos(sqrt(1-z(i)^2)*t)),'color',.8*colors(i,:),'linewidth',1)
    elseif z(i) == 1
        sig = -z(i);
        A = x0;
        B = -x0*sig;
        plot(t,exp(sig*t).*(A*t + B),'color',.8*colors(i,:),'linewidth',1)
    else
        sig1 = -(z(i) + sqrt(z(i)^2-1));
        sig2 = -(z(i) - sqrt(z(i)^2-1));
        
        C = [sig1, sig2;1,1]^-1*[0;x0];
        A = C(1);
        B = C(2);
        plot(t,A*exp(sig1*t) + B*exp(sig2*t),'color',.8*colors(i,:),'linewidth',1)
    end
    
    subplot(spX)
    hold on
    plot(r,1./sqrt((1 - r.^2).^2 + (2*z(i).*r).^2),'color',.8*colors(i,:),'linewidth',1)
    
    subplot(spF)
    hold on
    plot(r,sqrt((1 - r.^2).^2 + (2*z(i).*r).^2),'color',.8*colors(i,:),'linewidth',1)
    
    subplot(spPhi)
    hold on
    plot(r,atan2d(2*z(i).*r,1 - r.^2),'color',.8*colors(i,:),'linewidth',1)
    
end


legCell = cell(length(z),1);
for i=1:length(z)
%     legCell{i} = sprintf('\zeta = %3.3f',z(i));
    legCell{i} = ['\zeta = ',num2str(z(i))];
end

subplot(spXofT)
hold on
legend(legCell,'location','southeast')

h.Position(3) = 650;