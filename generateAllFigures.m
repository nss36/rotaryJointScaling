
%Nicholas Szczecinski
%West Virginia University
%13 October 2020
    
clearvars
close all

tPowerMin = -3;
tPowerMax = 2;
lPowerMin = -4;
lPowerMax = 1;

m0 = 12;
c0 = 1.31e3; %2.853e3; %see zetaAcrossScales.m
k0 = 12e3; %100e3;
s = sqrt(1e-3);

baseParams = [m0,c0,k0,s];

Tn = 2*pi*sqrt(m0/k0);

%List of labels, characteristic times, and characteristic lengths.
%Name, period of locomotion (or other behavior), limb length, t-axis label placement, L-axis label placement.
            
behaviorCell = {'Drosophila',       [50e-3, 130e-3],    0.5e-3, 1.5,    1;... %Wosnitza et al. 2013
                'Mus',              [.1, .33],          2e-2,   1.5,    1;... %Herbin et al. 2007
                'Rattus',           [.2, .6],           5e-2,   1.5,    1;... %Hruska et al. 1979
                'Aplysia grasp',    [1.25, 5],          1.5e-2, 1.5,    1;... %Shaw et al. 2015
                'Struthio'          [.25, 1],           1.26,   1.5,    1;... %Daley, Channon, et al. 2016
                'Homo finger tap',  [35e-3, 1],         6e-2,   1.5,    1;... %Aoki et al. 2003, Aoki et al. 2005
                'Felix',            [250e-3, 700e-3],   25e-2,  1.5,    1;... %Grillner 1975
                'Homo',             [500e-3, 1],        1,      1.5,    1;... %Grillner et al. 1979
                'Carausius',        [.6, 1.8],          4e-2,   1.5,    1;... %Cruse and Bartling 1995
                'Blaberus',         [.067, .25],        1.3e-2, 1.5,    1;... %Bender et al. 2010
                'Equus',            [.44, .9],          1.5,    1.5,    1;... %Hildebrand 1959, Hooper et al. 2009
                'Periplaneta',      [.0435, 0.667],     8e-3,   1.5,    1;... %Delcomyn 1971
                'Cicindela',        [.0556, .2],        .33e-2, 1.5,    1;... %Haselsteiner et al. 2014
                'Elephantidae',     [.5, 3],            1.7,    1.5,    1;... %Hutchinson et al. 2006 
%                 'Paratarsotomus',   [.0111, .1],        1e-3,   1.5,    1}; %Wu, Wright et al. 2010, 50 degrees Celsius 
                'Paratarsotomus',   [.05, .1],        1e-3,   1.5,    1;... %Wu, Wright et al. 2010, 30 degrees Celsius  
                'Osphranter',       [1/2.2,1/2.2+.1],      0.58,   1.5,    1;...
                'Petrogale',        [1/2.8,1/2.8+.1],      0.33,   1.5,    1;...
                'Pedetes',          [1/3,1/3+.1],          0.1,    1.5,    1;...
                'Potoroidae',       [1/3.5,1/3.5+.1],      0.08,   1.5,    1;...
                'Notomys',          [1/6,1/6+.1],          0.03,   1.5,    1};
  
behaviorCell = {'fruit fly',            [50e-3, 130e-3],    0.5e-3, 1.5,    1;... %Wosnitza et al. 2013
                'mouse',                [.1, .33],          2e-2,   1.5,    1;... %Herbin et al. 2007
                'rat',                  [.2, .6],           5e-2,   1.5,    1;... %Hruska et al. 1979
                'sea slug (grasp)',     [1.25, 5],          1.5e-2, 1.5,    1;... %Shaw et al. 2015
                'emu'                   [.25, 1],           1.26,   1.5,    1;... %Daley, Channon, et al. 2016
                'human (finger tap)',   [35e-3, 1],         6e-2,   1.5,    1;... %Aoki et al. 2003, Aoki et al. 2005
                'cat',                  [250e-3, 700e-3],   25e-2,  1.5,    1;... %Grillner 1975
                'human',                [500e-3, 1],        1,      1.5,    1;... %Grillner et al. 1979
                'stick insect',         [.6, 1.8],          4e-2,   1.5,    1;... %Cruse and Bartling 1995
                'discoid cockroach',    [.067, .25],        1.3e-2, 1.5,    1;... %Bender, Simpson et al. 2011
                'horse',                [.44, .9],          1.5,    1.5,    1;... %Hildebrand 1959, Hooper et al. 2009
                'American cockroach',   [.0435, 0.667],     8e-3,   1.5,    1;... %Delcomyn 1971, The Locomotion of the Cockroach Periplaneta Americana
                'tiger beetle',         [.0556, .2],        .33e-2, 1.5,    1;... %Haselsteiner et al. 2014
                'elephant',             [.5, 5],            1.7,    1.5,    1;... %Hutchinson et al. 2006 
                'mite',                 [.025, .2],        1e-3,   1.5,    1;... %Wu, Wright et al. 2010, 30 degrees Celsius  
                'red kangaroo',   1/2.2*[.95, 1.05],      0.58,   1.5,    1;...
                'spring hare',      1/3*[.95, 1.05],          0.1,    1.5,    1;...
                'rat kangaroo',   1/3.5*[.95, 1.05],      0.08,   1.5,    1;...
                'hopping mouse',    1/6*[.95, 1.05],          0.03,   1.5,    1};

            
nSampsPlot = 200;
nSampsInvDyn = 5;

%Plot the scaling figure used in the NeuroNex proposal. The understanding
%behind it is incomplete, but it represents an important development in
%this project.
hNeuroNex = NEURONEXscalingPlot(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorCell,baseParams,nSampsPlot);

%Plot how the forces that contribute to the steady-state response scale
%with L and T. This plot is a lot to take in, but it gives a useful
%summary.
hSteadyStateResponse = scalingInverseDynamicsAcrossScales(tPowerMin,tPowerMax,lPowerMin,lPowerMax,baseParams,nSampsInvDyn);

%Plot how the phase lag of the limb displacement, phi, scales with L and T.
%This is a proxy for where the actuator energy goes during steady-state
%motion.
%Plot how the damping ratio zeta scales with L and T. This is a proxy for 
%where unwanted energy from perturbations goes.
%Plot the size and cycle duration of various animal behaviors on the L and
%T plot. This shows that animal behaviors span several dynamic regimes, and
%some animals operate within several.
[g1,g2] = phiAndX(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorCell,baseParams,nSampsPlot);

keyboard
%Show that the phase angle predicts the full continuum of observed motor
%patterns that control reaching motions.
t = 0:.01:5;
x = .5+.5*erf(sqrt(pi) * (t - 2.5));
k = phiGeneratesResponse([0,63,90,117,180],x,t);

%Plot a color key to go with this figure.
colorVec = autumn(181);
figure;
for i=1:181
    polarplot((i-1)/180*pi*[1,1],[0,1],'color',colorVec(i,:))
    hold on
end
pax = gca;
pax.RTick = [];
pax.ThetaLim = [0,180];

%Plot system responses, mostly just for fun. The simulateJointResponse
%function is useful for this.
%QS, overdamped
L = 1e-2;
T = 1;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Quasi-static, overdamped')
%QS, underdamped
L = 1.7;
T = 10;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Quasi-static, underdamped')
%Inertial, underdamped
L = 1.7;
T = 0.5;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Inertial, underdamped')
%viscous, overdamped
L = 1e-2;
T = 1e-2;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Viscous, overdamped')
%Inertial, overdamped
L = .03;
T = 10^-2.8;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Inertial, overdamped')
