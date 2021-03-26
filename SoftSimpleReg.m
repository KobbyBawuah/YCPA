winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight


%physics constants 
dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15%time to sim
f = 230e12;%freq
lambda = c_c/f;%wavelength

xMax{1} = 20e-6; %setting up region
%steps
nx{1} = 200;    
ny{1} = 0.75*nx{1};


Reg.n = 1;

%creating the mu matrix
mu{1} = ones(nx{1},ny{1})*c_mu_0;

% Create inclusion
epi{1} = ones(nx{1},ny{1})*c_eps_0;

%epi{1}(125:150,55:95)= c_eps_0*11.3;
%first inclusion
epi{1}(round(125*nx{1}/200):round(150*nx{1}/200),...
    round(55*nx{1}/200): round(95*nx{1}/200))= c_eps_0*11.3
%second inclusions structure
%comment above and unccomment below to change structure
% epi{1}(125:150,55:95)= c_eps_0*11.3;
% epi{1}(125:150,1:41) = c_eps_0*11.3;
% epi{1}(125:150,110:150) = c_eps_0*11.3;
% 
% epi{1}(50:60,20:25) = c_eps_0*11.3;
% epi{1}(140:150,20:25) = c_eps_0*11.3;
%     
% epi{1}(50:60,125:130) = c_eps_0*11.3;
% epi{1}(125:135,125:130) = c_eps_0*11.3;
% epi{1}(140:150,125:130) = c_eps_0*11.3;
%     
% epi{1}(140:150,30:35) = c_eps_0*11.3;
% epi{1}(140:150,40:45) = c_eps_0*11.3;
% epi{1}(140:150,110:115) = c_eps_0*11.3;

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1}; %Step size
dt = 0.25*dx/c_c; %Time step
nSteps = round(tSim/dt*2); %Simulation length
yMax = ny{1}*dx; %Y dimension size
nsteps_lamda = lambda/dx %Number of lambdas

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];


bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
%wave param
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = -0.05; %st = 15e-15;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

bc{1}.s(2).xpos = nx{1}/(4) + 1;
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
%Wave parameters
mag = 5;
phi = pi;
omega = f*2*pi;
betap = 0;
t0 = 45e-15;
%st = 15e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 2*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; %Boundary conditions

Plot.y0 = round(y0/dx);

% The type of boundary
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

bc{2}.xm.type = 'a';
bc{2}.xp.type = 'a';
bc{2}.ym.type = 'a';
bc{2}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg




