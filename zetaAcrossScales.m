close all

%Nicholas Szczecinski
%West Virginia University
%13 October 2020

%Values from the literature:
%There are six studies here, in order:
%Garcia, Kuo, Peattie, Wang, and Full 2000, cockroach tibia
%Zakotnik, Matheson, and Duerr 2006, locust tibia
%Hajian and Howe 1997, human finger
%Weiss, Hunter, and Kearney 1988, human ankle
%Stein, Zehr, Lebiedowska, Popovic, Sheiner, and Chizeck 1996, human tibia
%Stein, Zehr, Lebiedowska, Popovic, Sheiner, and Chizeck 1996, human leg
%(estimated by Szczecinski based on the logarithmic decrement in the data
%in Figure 6, control subject).
L =         [.0113, .0228,  .0976,     .22,      0.50,      .91];
zeta =      [6.5,	3.14,   1,          1/3,    0.098,      .05];
zetaUp =    [6.61,  9,      1.41,      .4,      0.219,      .06];
zetaDn =    [6.39,  .3,     .43,      .25,      0.041,      .04];

k0 = 12e3; %100e3;
m0 = 12;
s = sqrt(1e-3);
g = 10;

lw = 1; %Linewidth in plots

%Just plot the "mean" values of zeta from each study against L to see how
%they scale. We expect roughly zeta \propto L^-1
h1 = figure;
plot(L,zeta,'linewidth',lw)
hold on
title('Average \zeta measurements across scales')
grid on
xlabel('L (m)')
ylabel('\zeta')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';


h2 = figure;
errorbar(L,zeta,zeta-zetaDn,zeta-zetaUp,'o','linewidth',lw)
hold on
title('Damping ratio follows rough \zeta \propto L^{-1} relationship')
grid on
xlabel('L (m)')
ylabel('\zeta')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlim([1e-2,1e0])
ylim([1e-2,1e1])
ax.Position = [0.2750    0.1106    0.4500    0.8000];

h3 = figure;
copyobj(ax,h3)

%Use nonlinear optimization to pick the value of c0 that best fits the
%data, based on k0, m0, and scaling laws. In this case, we assume that only
%elastic potential energy is stored in the limb, not gravitational
%potential energy.
fElasOnly = @(c0) sum(( log(zeta) - log(c0*s^2/(2*sqrt(k0*s^2*1/3*m0))*1./L)).^2);

lb = 0;
ub = inf;
c0elasOnly = fmincon( fElasOnly,.1,[],[],[],[],lb,ub);

zetaElasOnly = c0elasOnly*s^2/(2*sqrt(k0*s^2*1/3*m0));

figure(h2)
hold on
pp1 = plot(L,zetaElasOnly./L,'linewidth',lw);
legend(pp1,'no grav potential energy')

%But what if k also depends on gravity?
%Before tuning this relationship with the animal data, let's just compute
%what we would expect zeta to be if k depended only on elastic elements,
%if it depended only on gravity, or if it depended on both.
l = logspace(-3,1);
kElas = k0*s^2*l.^3;
kGrav = m0*g/2*l.^4;
kTotal = kElas + kGrav;
J = 1/3*m0*l.^5;
c = 2*c0elasOnly*s^2*l.^3;

zetaGrav = c./(2*sqrt(kGrav.*J));
zetaElas = c./(2*sqrt(kElas.*J));
zetaTotal = c./(2*sqrt(kTotal.*J));

figure
plot(l,zetaGrav,'linewidth',lw)
hold on
plot(l,zetaElas,'linewidth',lw)
plot(l,zetaTotal,'linewidth',lw)
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
legend('k = k_{grav}','k = k_{r,elas}','k = k_{r,grav} + k_{r,elas}')
ylabel('\zeta')
xlabel('L (m)')
title('Comparing how \zeta scales with different forms of k')

%Now, let us use the same nonlinear least squares process as before, but
%this time include gravitational potential energy in our formulation.
fzeta = @(c0,L) c0.*s^2.*L.^3./(2.*sqrt( (k0.*s^2.*L.^3 + m0.*10/2.*L.^4).*1/3.*m0.*L.^5   ));
f3 = @(c0) sum(( log(zeta) - log(fzeta(c0,L))).^2 );
lb = 0;
ub = inf;
c0 = fmincon( f3,.1,[],[],[],[],lb,ub);

%Add the zeta curve to the figure of empirical data points.
figure(h3)
hold on
plot(l,fzeta(c0,l),'linewidth',lw)

%How do elastic and gravitational energy trade off as a function of scale,
%anyway? Plot each of the individual stiffness terms along with their sum,
%to observe how they differ over scale.
h4 = figure;
plot(l,kGrav,'linewidth',lw)
hold on
plot(l,kElas,'linewidth',lw)
plot(l,kTotal,'linewidth',lw)
legend('k_{r,grav}','k_{r,elas}','k_{r,grav} + k_{r,elas}','location','southeast')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('L (m)')
ylabel('stiffness (Nm/rad)')
title('The balance between elastic and gravitational potential depends on scale')

%How do these different terms for k alter the natural period? Plot this
%dividing line as we would see it on the T and L plot, in all 3 cases for
%comparison.
colorVec = lines(6);
h5 = figure;
plot(l,2*pi./sqrt(kGrav./J),'color',colorVec(6,:),'linewidth',lw)
hold on
plot(l,2*pi./sqrt(kElas./J),'color',colorVec(2,:),'linewidth',lw)
plot(l,2*pi./sqrt(kTotal./J),'color',colorVec(4,:),'linewidth',lw)
legend('T_n when k = k_{r,grav}','T_n when k = k_{r,elas}','T_n when k = k_{r,grav} + k_{r,elas}','location','southeast')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('L (m)')
ylabel('T (s)')
title('The natural period depends on which terms are included in potential energy')
ylim([1e-3,1e1])

%What is the relative contribution of k_elas and k_grav to k?
h6 = figure;
plot(l,kGrav./kTotal,'linewidth',lw)
hold on
plot(l,kElas./kTotal,'linewidth',lw)
ax = gca;
ax.XScale = 'log';
xlabel('L (m)')
ylabel('')
legend('k_{r,grav}/k_{r,total}','k_{r,elas}/k_{r,total}','location','east')
title('The joint stiffness of larger animals is due to gravity')
xlim([1e-3,1e1])

fprintf('The baseline viscous damping coefficient c0 = %3.3e.\n',c0)
