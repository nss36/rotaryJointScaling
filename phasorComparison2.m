m = [2,2,10,100]*.1;
c = [2,10,2,20]*.1;
k = [10,2,2,20]*.1;

m = [2,2,10,2,2]*.1;
c = [2,10,2,10,10]*.1;
k = [10,2,2,2,2]*.1;

omega_n = sqrt(k./m);
omega = [omega_n(2),omega_n(2),omega_n(2),.05*omega_n(2),20*omega_n(2)];

nrows = 5;
ncols = 3;

X = 1./sqrt( (k - m.*omega.^2).^2 + (c.*omega).^2 );

phi = atan2(c.*omega,k - m.*omega.^2);

T = 2*pi./omega;


h = figure;

clines = lines(7);
cmap = [clines(5,:);clines(7,:);clines(6,:)];
alpha = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'};

for i=1:nrows
    %Bar graph of parameter values
    subplot(nrows,ncols,ncols*(i-1)+1)
%     bar(1:3,[m(i),c(i),k(i)])
    bar(1,m(i),'FaceColor',cmap(3,:),'EdgeAlpha',0)
    hold on
    bar(2,c(i),'FaceColor',cmap(2,:),'EdgeAlpha',0)
    bar(3,k(i),'FaceColor',cmap(1,:),'EdgeAlpha',0)
    xticks([1,2,3])
    xticklabels({'m','c','k'})
    title([alpha{ncols*(i-1)+1}','. Parameter values'])
    
    %Vector representation
    subplot(nrows,ncols,ncols*(i-1)+2)
    xmax = 1.25;
    
    thsamp = linspace(0,pi);
    plot(cos(thsamp),sin(thsamp),'-.','color',[.5,.5,.5]);
    hold on    
    arrF = plot([0,cos(phi(i))],[0,sin(phi(i))],'k--');
    arrK = plot([0,k(i)*X(i)],[0,0],'color',cmap(1,:));
    arrC = plot([k(i)*X(i),k(i)*X(i)],[0,c(i)*omega(i)*X(i)],'color',cmap(2,:));
    arrM = plot([k(i)*X(i),k(i)*X(i)-m(i)*omega(i)^2*X(i)],[c(i)*omega(i)*X(i),c(i)*omega(i)*X(i)],'color',cmap(3,:));    
    axis equal
    xlim([-xmax,xmax])
    ylim([0,xmax])
    title([alpha{ncols*(i-1)+2}','. Phasor plot, r = ',num2str(omega(i)/omega_n(i))])
    if i == nrows
        legend('','applied force','elas/grav force','viscous force','inertial force','location','southeast')
    end
    
    %Response
    t = linspace(0,3*T(i),300);
    subplot(nrows,ncols,ncols*(i-1)+3)
    plot(t,cos(omega(i)*t),'k--')
    hold on
    plot(t,X(i)*cos(omega(i)*t-phi(i)),'color',cmap(1,:))
    title([alpha{ncols*(i-1)+3}','. Steady state response'])
    if i == nrows
        legend('applied force','displacement','location','southeast')
    end
    xlabel('time (s)')
    ylabel('force or displacement')
end

h.Position(2) = 0;
h.Position(3) = 650;
h.Position(4) = 750;
