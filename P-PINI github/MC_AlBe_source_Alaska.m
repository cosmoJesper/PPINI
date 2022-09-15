% function MC_AlBe_source

% Source sample generator for the P-PINI code.
%
% The code will generate a .mat file containing:
%   10Be concentrations [atoms/gr]
%   26Al concentrations [atoms/gr]
%   All input parameters
%
% The workflow of the program is:
%
% ->      Run MC_AlBe_source.m
%                   |
%          Run MC_AlBe_sink.m
%                   |
%         Run PPINI_analylsis.m

% Name the output source file and determine the Number of simulations you want:
name = 'mc_source_Alaska.mat';                               % Filename
Nsim = 1e7;                                                  % Number of simulations

%*********** Geographical settings ************
lat = 63.9;                                                  % Latitude of sink site
lon = -148.9;                                                % Longitude of sink site, negative if W
sink_elev = 563;                                             % Elevation of sampling site [m asl.]
max_elev = 5000;                                             % Maximum elevation of the possible source area [m asl.]

%*********** Landscape settings **************
rho_min = 2.5;                                               % Minimum rock density in source [g/cm3]
rho_max = 2.7;                                               % Maximum rock density in source [g/cm3]
min_ero = 10;                                                % Minimum source interglacial erosion rate [m/Myr]
max_ero = 1500;                                              % Miximum source interglacial erosion rate [m/Myr]

%************* Advanced settings ***************
P_SLHL_Be = 4;                                                  % 10Be production sea-level high-latitude [at/g/yr]
R = 6.97;                                                    % 26Al/10Be production ratio at sea-level high-latitude
P_SLHL_Al = R*P_SLHL_Be;
egla = 0e-6;                                                 % Glacial erosion rate [m/yr]
pluck_min = 0.01;                                            % Minimum plucking depth pr. glaciation [m] 
pluck_max = 10;                                              % Maximum plucking depth pr. glaciation [m] 
ice_lengths = 10000;                                         % Maximum potential length of a glaction at site
n_ice_max = 3;                                               % Number of glaciations possible

dt = 2000;                                                   % Timesteps in calculations [yr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           Code form here           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete(gcp('nocreate'))
cd other/code
close all
rng(1)

nn = 100000;

n_elev = 30;
elev = linspace(sink_elev,max_elev,n_elev);                       % Elevation of possible source area [m asl.]
eint_mc = logspace(log10(min_ero*1e-6),log10(max_ero*1e-6),nn);   % Erosion rates of source area [m/Myr]
rho_pool = linspace(rho_min,rho_max,nn);                                  % Rock density in source [g/cm3]
pluck_depth = logspace(log10(pluck_min),log10(pluck_max),nn);   % Plucking depth pr. glaciation [m]

%************* Setup calculations **************
elev_rough = linspace(min(elev),max(elev),n_elev);
Pspal_Be_rough = zeros(1,n_elev);
Pspal_Al_rough = zeros(1,n_elev);

%*********** Cosmo constants ************

Lspal_Be = 160;                                         % 10Be spallation attenuation length [g/cm^2]
Lspal_Al = 160;                                         % 26Al spallation attenuation length [g/cm^2]
Lnmc_Be = 1500;                                         % 10Be neutron-capture attenuation length [g/cm^2]
Lnmc_Al = 1500;                                         % 26Al neutron-capture attenuation length [g/cm^2]
Lfm_Be = 4320;                                          % 10Be fast muon attenuation length [g/cm^2]
Lfm_Al = 4320;                                          % 26Al fast muon attenuation length [g/cm^2]

TBe = 1.387e6;                                          % 10Be half-life [yr]
TAl = 0.705e6;                                          % 26Al half-life [yr]
lambda_Be = log(2)/TBe;                                 % 10Be mean lifetime [yr^-1]
lambda_Al = log(2)/TAl;                                 % 26Al mean lifetime [yr^-1]

%%%
gmr = -0.03417; % Assorted constants
dtdz = 0.0065; % Lapse rate from standard atmosphere
load consts_LSD
z = linspace(0,max_ero*1e-6*50*dt,10000);
z_rho_max = z*100*max(rho_pool);

f10C = 0.704;
f10D = 0.1828;
c10.fC = f10C;
c10.fD = f10D;
c10.Natoms = consts.Natoms10;
c10.fstar = 0.00157;
c10.sigma190 = 37.8e-30;

f26C = 0.296;
f26D = 0.6559;
c26.fC = f26C;
c26.fD = f26D;
c26.Natoms = consts.Natoms26;
c26.fstar = 0.0118;
c26.sigma190 = 521e-30;

fitfun = fittype( @(a1,e1,a2,e2,a3,e3,x) a1*exp(x*e1) + a2*exp(x*e2) + a3*exp(x*e3));
x0_Be = [1 -1/Lnmc_Be 0.1 -1/Lfm_Be 0.01 -1/Lfm_Be]; 
x0_Al = [1 -1/Lnmc_Al 0.1 -1/Lfm_Al 0.01 -1/Lfm_Al]; 

exp_vars = zeros(12,n_elev);

options_Be = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0 -10 0 -10 0 -10], 'Upper',[10 0 10 0 10 0],'StartPoint',x0_Be);
options_Al = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0 -10 0 -10 0 -10], 'Upper',[10 0 10 0 10 0],'StartPoint',x0_Al); 

global h p
h = waitbar(0, 'Calculating production rates, please wait...');
p = 1;                 

tic

for i = 1:n_elev
    
    LSD_DG(lat, lon, elev_rough(i),0,1e7,-1,10);               % Calculating 10Be production scaling
    load LSDout;                                               % Obtaining 10Be production scaling
    P_time = P_SLHL_Be*mean(LSDout.Be);                           % Calculating 10Be production at height
    Pspal_Be_rough(i) = P_time(1);                             
    
    LSD_DG(lat, lon, elev_rough(i),0,1e7,-1,26);            % Calculating 26Al productions scaling
    load LSDout;                                            % Obtaining 26Al production scaling
    P_time = P_SLHL_Al*mean(LSDout.Al);                      % Calculating 26Al production at height
    Pspal_Al_rough(i) = P_time(1);
    
    pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev(i).*dtdz)) ) );
    
    P_mu_Be = P_mu_total_edited(z_rho_max,pressure,c10,'no');
    [fitted_curve_Be,gof] = fit(z_rho_max',P_mu_Be',fitfun,options_Be);
    exp_vars(1,i) = fitted_curve_Be.a1;
    exp_vars(2,i) = fitted_curve_Be.e1;
    exp_vars(3,i) = fitted_curve_Be.a2;
    exp_vars(4,i) = fitted_curve_Be.e2;
    exp_vars(5,i) = fitted_curve_Be.a3;
    exp_vars(6,i) = fitted_curve_Be.e3;
    
    P_mu_Al = P_mu_total_edited(z_rho_max,pressure,c26,'no');
    [fitted_curve_Al,gof] = fit(z_rho_max',P_mu_Al',fitfun,options_Al);
    exp_vars(7,i) = fitted_curve_Al.a1;
    exp_vars(8,i) = fitted_curve_Al.e1;
    exp_vars(9,i) = fitted_curve_Al.a2;
    exp_vars(10,i) = fitted_curve_Al.e2;
    exp_vars(11,i) = fitted_curve_Al.a3;
    exp_vars(12,i) = fitted_curve_Al.e3;
      
     t = toc;
     waitbar(i/n_elev, h, ['Calculating production rates, please wait...: ',num2str(round(t/(i/n_elev)-t)),' seconds']);
    
end

close(h)

Pspal_Be_all = interp1(elev_rough,Pspal_Be_rough,elev);     % 10Be production
Pspal_Al_all = interp1(elev_rough,Pspal_Al_rough,elev);     % 26Al production

prob = 1/nn;                                                % Probability
cumprob = prob:prob:1;                                      % Cummulative probability for elevations

source_elev = elev(elev>=sink_elev);
source_prob = cumprob(elev>=sink_elev);
source_prob = source_prob/max(source_prob);                 
Pspal_Be_all = interp1(elev,Pspal_Be_all,source_elev);      % 10Be production by spallation [atoms/(g*yr)]
Pnmc_Be_all = Pspal_Be_all*0.015;                           % 10Be production by neutron-capture [atoms/(g*yr)]
Pfm_Be_all = Pspal_Be_all*0.005;                            % 10Be production by fast muons [atoms/(g*yr)]   
Pspal_Al_all = interp1(elev,Pspal_Al_all,source_elev);      % 26Al production by spallation [atoms/(g*yr)]
Pnmc_Al_all = Pspal_Al_all*0.018;                           % 26Al production by neutron-capture [atoms/(g*yr)]
Pfm_Al_all = Pspal_Al_all*0.006;                            % 26Al production by fast muons [atoms/(g*yr)]

% Matrices for saving at the end
NBe_fin = zeros(1,Nsim);                        % Final concentrations 10Be pre-sink burial
NAl_fin = zeros(1,Nsim);                        % Final concentrations 26Al pre-sink burial
Tice = zeros(1,Nsim);                           % Ice covered time           
source_pluck = zeros(1,Nsim);                   % Source total plucked
Pspal = zeros(1,Nsim);                          % 10Be production rate
NBe_source = zeros(1,Nsim);                     % 10Be inherited from source
NAl_source = zeros(1,Nsim);                     % 26Al inherited from source
erate_source = zeros(1,Nsim);                   % Erosion rates in source
elevation = zeros(1,Nsim);                      % Source elevation
max_burial = zeros(1,Nsim);

%%% Creating a waitbar
global h p
D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);
p = 1;

parpool;            % Start parallel workers
tic                 % Start timer

parfor (j=1:Nsim,3)

    rho = rho_pool(randperm(nn,1));

    % Defining ice cover steps (each step represents dt=2000 yr)
    % pfac = 0 means ice covered, pfac = 1 means ice free
    ice = randi(round(ice_lengths/dt),n_ice_max*2,1)-1;
    ticker = 1;
    pfac = zeros(1,sum(ice));
    
    for i=1:2:n_ice_max*2
        ticker = ticker + ice(i);
        pfac(ticker:ticker+ice(i+1)) = 1;
        ticker = ticker + ice(i+1);
    end
    
    pfac = [pfac ones(1,50)];
    ntime = length(pfac);
    
    I_ice = find(pfac == 0);                       
    Tice(j) = dt*1e-6*length(I_ice);                % Total time of ice cover [Myr]    
    
    eint = eint_mc(randperm(nn,1));                 % Picking a random erosion rate

    %%% Following a elevation profile
    r = rand;                                       % Pick random productrion rate
    idx = find(r < source_prob);                    
    prod = idx(1);                               
    elevation(j) = source_elev(prod);
    Pspal_Be = Pspal_Be_all(prod);
    Pnmc_Be = Pnmc_Be_all(prod);
    Pfm_Be = Pfm_Be_all(prod);
    Pspal_Al = Pspal_Al_all(prod);
    Pnmc_Al = Pnmc_Al_all(prod);
    Pfm_Al = Pfm_Al_all(prod);

    a1_Be = exp_vars(1,prod);
    e1_Be = exp_vars(2,prod);
    a2_Be = exp_vars(3,prod);
    e2_Be = exp_vars(4,prod);
    a3_Be = exp_vars(5,prod);
    e3_Be = exp_vars(6,prod);
    
    a1_Al = exp_vars(7,prod);
    e1_Al = exp_vars(8,prod);
    a2_Al = exp_vars(9,prod);
    e2_Al = exp_vars(10,prod);
    a3_Al = exp_vars(11,prod);
    e3_Al = exp_vars(12,prod);

    burial = zeros(ntime,1);                        % List for saving burial depth of sample
    burial(1) = 0;                                  % Final burial depth

    for i=2:ntime-1                                             % Looping to calculate burial depth with time (stepping backwards in time)
        if pfac(i) == 0                                         % If ice cover
            if pfac(i-1) == 0 && pfac(i+1) == 1                 % If start of ice age 
                pluck_now = pluck_depth(randperm(nn,1));        % Random plucking amount
                burial(i) = burial(i-1) + pluck_now;            % Do random plucking
                source_pluck(j) = source_pluck(j) + pluck_now;  % Add to total plucking amount
            else
                burial(i) = burial(i-1) + dt*egla; 
            end
        else                                                    % If not ice covered
            burial(i) = burial(i-1) + dt*eint;                  % Use interglacial erosion rate
        end
    end
    burial(ntime) = burial(ntime-1) + dt*eint;                  % Last time step
    
    % Calculating steady state starting point
%     NBe = zeros(1,ntime);
%     NAl = zeros(1,ntime);
%     fBe = lambda_Be + rho*eint*100/Lspal_Be;
%     fAl = lambda_Al + rho*eint*100/Lspal_Al;
%     NBe(ntime) = Pspal_Be*exp(-rho*100*burial(ntime)/Lspal_Be)/fBe;
%     NAl(ntime) = Pspal_Al*exp(-rho*100*burial(ntime)/Lspal_Al)/fAl;
%     fBe = lambda_Be + rho*eint*100/Lnmc_Be;
%     fAl = lambda_Al + rho*eint*100/Lnmc_Al;
%     NBe(ntime) = NBe(ntime) + Pnmc_Be*exp(-rho*100*burial(ntime)/Lnmc_Be)/fBe;
%     NAl(ntime) = NAl(ntime) + Pnmc_Al*exp(-rho*100*burial(ntime)/Lnmc_Al)/fAl; 
%     fBe = lambda_Be + rho*eint*100/Lfm_Be;
%     fAl = lambda_Al + rho*eint*100/Lfm_Al;
%     NBe(ntime) = NBe(ntime) + Pfm_Be*exp(-rho*100*burial(ntime)/Lfm_Be)/fBe;
%     NAl(ntime) = NAl(ntime) + Pfm_Al*exp(-rho*100*burial(ntime)/Lfm_Al)/fAl;
%     
    %%%%%%%%%%    
    NBe = zeros(1,ntime);
    NAl = zeros(1,ntime);
    fBe = lambda_Be + rho*eint*100/Lspal_Be;
    fAl = lambda_Al + rho*eint*100/Lspal_Al;
    NBe(ntime) = Pspal_Be*exp(-rho*100*burial(ntime)/Lspal_Be)/fBe;
    NAl(ntime) = Pspal_Al*exp(-rho*100*burial(ntime)/Lspal_Al)/fAl;
    fBe = lambda_Be - rho*eint*100*e1_Be;
    fAl = lambda_Al - rho*eint*100*e1_Al;
    NBe(ntime) = NBe(ntime) + a1_Be*exp(rho*100*burial(ntime)*e1_Be)/fBe;
    NAl(ntime) = NAl(ntime) + a1_Al*exp(rho*100*burial(ntime)*e1_Al)/fAl; 
    fBe = lambda_Be - rho*eint*100*e2_Be;
    fAl = lambda_Al - rho*eint*100*e2_Al;
    NBe(ntime) = NBe(ntime) + a2_Be*exp(rho*100*burial(ntime)*e2_Be)/fBe;
    NAl(ntime) = NAl(ntime) + a2_Al*exp(rho*100*burial(ntime)*e2_Al)/fAl;
    fBe = lambda_Be - rho*eint*100*e3_Be;
    fAl = lambda_Al - rho*eint*100*e3_Al;
    NBe(ntime) = NBe(ntime) + a3_Be*exp(rho*100*burial(ntime)*e3_Be)/fBe;
    NAl(ntime) = NAl(ntime) + a3_Al*exp(rho*100*burial(ntime)*e3_Al)/fAl;
    %%%%%%
    
    
    for kk=(ntime-1):-1:1 % Go through all the burial timesteps from oldset to now.

          Pz_Be = pfac(kk)*(Pspal_Be*exp(-rho*100*burial(kk)/Lspal_Be)...           % Balco
                +a1_Be*exp(e1_Be*rho*100*burial(kk))...
                +a2_Be*exp(e2_Be*rho*100*burial(kk))...
                +a3_Be*exp(e3_Be*rho*100*burial(kk)));

          Pz_Al = pfac(kk)*(Pspal_Al*exp(-rho*100*burial(kk)/Lspal_Al)...           % Balco
                +a1_Al*exp(e1_Al*rho*100*burial(kk))...
                +a2_Al*exp(e2_Al*rho*100*burial(kk))...
                +a3_Al*exp(e3_Al*rho*100*burial(kk)));


%             Pz_Be = pfac(kk)*(Pspal_Be*exp(-rho*100*burial(kk)/Lspal_Be)...       %Heisinger
%                 +Pnmc_Be*exp(-rho*100*burial(kk)/Lnmc_Be)...
%                 +Pfm_Be*exp(-rho*100*burial(kk)/Lfm_Be));
% 
%             Pz_Al = pfac(kk)*(Pspal_Al*exp(-rho*100*burial(kk)/Lspal_Al)...         %Heisinger
%                 +Pnmc_Al*exp(-rho*100*burial(kk)/Lnmc_Al)...
%                 +Pfm_Al*exp(-rho*100*burial(kk)/Lfm_Al));

        % Calculate new concentrations with decays of nuclides
        NBe(kk) = NBe(kk+1)*exp(-dt*lambda_Be) + Pz_Be*(1-exp(-dt*lambda_Be))/lambda_Be;

        NAl(kk) = NAl(kk+1)*exp(-dt*lambda_Al) + Pz_Al*(1-exp(-dt*lambda_Al))/lambda_Al;
    end

    NBe_fin(j) = NBe(1); % Newest concentrations
    NAl_fin(j) = NAl(1);

    Pspal(j) = Pspal_Be;
    NBe_source(j) = NBe(1); 
    NAl_source(j) = NAl(1); 
    erate_source(j) = eint;
    max_burial(j) = burial(end);
    
    if j==(Nsim*0.001)*floor(j/(Nsim*0.001))
        send(D, j);
    end   
        
end

delete(gcp('nocreate'));

waitbar(1,h,'Saving...')    

%%%%%%%%%%%%%%%%%% Saving results %%%%%%%%%%%%%%%%%%%%%%%

mc_source.NBe_fin = NBe_fin;
mc_source.NAl_fin = NAl_fin;
mc_source.Nsim = Nsim;
mc_source.storage_time = zeros(size(NBe_fin));
mc_source.storage_depth = zeros(size(NBe_fin));
mc_source.source_pluck = source_pluck;
mc_source.Tice = Tice;
mc_source.P_source = Pspal;
mc_source.NBe_source = NBe_source;
mc_source.NAl_source = NAl_source;
mc_source.erate_source = erate_source;
mc_source.elevation_source = elevation;
mc_source.max_burial = max_burial;

save(name,'mc_source','-v7.3')

figure
plot(NBe_fin(1:1000:end),NAl_fin(1:1000:end),'.k')
hold on
title('Source simulation concentrations')
plot([min(NBe_fin) max(NBe_fin)],[min(NBe_fin) max(NBe_fin)]*R,'-r')
xlabel('^1^0Be [at/g]')
ylabel('^2^6Al [at/g]')

cd ..\..
close(h)
toc

% end


function nUpdateWaitbar(~)
    global p h
        t = toc;
        waitbar(p/1000, h, ['Running simulations, please wait... ',num2str(round(t/(p/1000)-t)),' seconds']);
        p = p + 1;
end

