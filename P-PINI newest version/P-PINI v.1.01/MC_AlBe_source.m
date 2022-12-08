function MC_AlBe_source

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
name = 'mc_source_Pulu.mat';                                 % Filename
Nsim = 5e6;                                                  % Number of simulations

%*********** Geographical settings ************
lat = 36;                                                  % Latitude of sink site
lon = 81.5;                                                % Longitude of sink site, negative if W
sink_elev = 2500;                                             % Elevation of sampling site [m asl.]
max_elev = 6000;                                             % Maximum elevation of the possible source area [m asl.]

%************* Landscape settings *************
rho_min = 2.6;                                               % Minimum rock density in source [g/cm3]
rho_max = 2.7;                                               % Maximum rock density in source [g/cm3]
min_ero = 10;                                                % Minimum source interglacial erosion rate [m/Myr]
max_ero = 1500;                                              % Miximum source interglacial erosion rate [m/Myr]

%************* Advanced settings **************
egla = 0e-6;                                                 % Glacial erosion rate [m/yr]
pluck_min = 0.01;                                            % Minimum plucking depth pr. glaciation [m] 
pluck_max = 10;                                              % Maximum plucking depth pr. glaciation [m] 
ice_lengths = 10000;                                         % Maximum potential length of a glaction at site [yr]
n_ice_max = 3;                                               % Number of glaciations possible (from 1 up to x)

P_SLHL_Be = 4;                                               % 10Be production sea-level high-latitude [at/g/yr]
R_SLHL = 6.8;                                                % 26Al/10Be production ratio at SLHL
P_SLHL_Al = R_SLHL*P_SLHL_Be;                                % 26Al production SLHL
dt = 2000;                                                   % Timesteps in calculations [yr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           Code form here           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete(gcp('nocreate')); cd other/code; close all; rng(1); tic

nn = 100000;
n_elev = 30;

%************* Setup for calculations **************
elev = linspace(sink_elev,max_elev,n_elev);                        % Elevation of possible source area [m asl.]
eint_mc = logspace(log10(min_ero*1e-6),log10(max_ero*1e-6),nn);    % Erosion rates of source area [m/Myr]
rho_pool = linspace(rho_min,rho_max,nn);                           % Rock density in source [g/cm3]
pluck_depth = logspace(log10(pluck_min),log10(pluck_max),nn);      % Plucking depth pr. glaciation [m]

elev_rough = linspace(min(elev),max(elev),n_elev);
Pspal_Be_rough = zeros(1,n_elev);
Pspal_Al_rough = zeros(1,n_elev);
z = linspace(0,min([1.9e5 max_ero*1e-6*50*dt*100*max(rho_pool)]),1000);

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
gmr = -0.03417;                                         % Assorted constants for atmospheric pressure
dtdz = 0.0065;                                          % Lapse rate from standard atmosphere
c10.fC = 0.704;                                         % 10Be specific constants
c10.fD = 0.1828;
c10.Natoms = 2.006e22;
c10.fstar = 0.00157;
c10.sigma190 = 37.8e-30;
c26.fC = 0.296;                                         % 26Al specific constants
c26.fD = 0.6559;
c26.Natoms = 1.003e22;
c26.fstar = 0.0118;
c26.sigma190 = 521e-30;

global h p
h = waitbar(0, 'Calculating production rates, please wait...');
p = 1;                 

%%%%%%%%%%%% Calcucation muon production rates %%%%%%%%%%%%%%
b_Be = zeros(6,n_elev);
b_Al = zeros(6,n_elev);
start_Be = [1 1 1 -1/Lnmc_Be -1/Lfm_Be -1/Lfm_Be];
start_Al = [1 1 1 -1/Lnmc_Al -1/Lfm_Al -1/Lfm_Al];

for i = 1:n_elev
    
    LSD_DG(lat, lon, elev_rough(i),0,1e7,-1,10);                % Calculating 10Be production scaling for spallation
    load LSDout;                                                % Obtaining 10Be production scaling for spallation
    P_time = P_SLHL_Be*mean(LSDout.Be);                         % Calculating 10Be production at height  for spallation
    Pspal_Be_rough(i) = P_time(1);                             
    
    LSD_DG(lat, lon, elev_rough(i),0,1e7,-1,26);                % Calculating 26Al productions scaling for spallation
    load LSDout;                                                % Obtaining 26Al production scaling for spallation
    P_time = P_SLHL_Al*mean(LSDout.Al);                         % Calculating 26Al production at height for spallation
    Pspal_Al_rough(i) = P_time(1);
    
    pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev(i).*dtdz)) ) ); % Atmospheric pressure at elevation
    
    P_mu_Be = P_mu_total_edited(z,pressure,c10,'no');           % Calculating 10Be muon production
    P_mu_Al = P_mu_total_edited(z,pressure,c26,'no');           % Calculating 26Al muon production
    
    objfcn = @(b,x) b(1)*exp(x*b(4)) + b(2)*exp(x*b(5)) + b(3)*exp(x*b(6)); % Multi-exponential function
    nlm_Be = fitnlm(z',P_mu_Be',objfcn,start_Be);               % Fitting multiexponential function to 10Be muon production
    b_Be(:,i) = nlm_Be.Coefficients.Estimate;
    start_Be = nlm_Be.Coefficients.Estimate; 
    nlm_Al = fitnlm(z',P_mu_Al',objfcn,start_Al);               % Fitting multiexponential function to 26Al muon production
    b_Al(:,i) = nlm_Al.Coefficients.Estimate;
    start_Al = nlm_Al.Coefficients.Estimate;
    
    t = toc;
    waitbar(i/n_elev, h, ['Calculating production rates, please wait...: ',num2str(round(t/(i/n_elev)-t)),' seconds']);
    
end

close(h)

Pspal_Be_all = interp1(elev_rough,Pspal_Be_rough,elev);     % 10Be spallation production
Pspal_Al_all = interp1(elev_rough,Pspal_Al_rough,elev);     % 26Al spallation production

prob = 1/n_elev;                                                % Probability
cumprob = prob:prob:1;                                      % Cummulative probability for elevations

% Matrices for saving at simulation information
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
    idx = find(r < cumprob);                    
    prod = idx(1);                               
    elevation(j) = elev(prod); 
    
    Pspal_Be = Pspal_Be_all(prod);
    Pspal_Al = Pspal_Al_all(prod);

    a1_Be = b_Be(1,prod);
    e1_Be = b_Be(4,prod);
    a2_Be = b_Be(2,prod);
    e2_Be = b_Be(5,prod);
    a3_Be = b_Be(3,prod);
    e3_Be = b_Be(6,prod);
    
    a1_Al = b_Al(1,prod);
    e1_Al = b_Al(4,prod);
    a2_Al = b_Al(2,prod);
    e2_Al = b_Al(5,prod);
    a3_Al = b_Al(3,prod);
    e3_Al = b_Al(6,prod);

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
    
    for kk=(ntime-1):-1:1 % Go through all the burial timesteps from oldset to now.

          Pz_Be = pfac(kk)*(Pspal_Be*exp(-rho*100*burial(kk)/Lspal_Be)...           % Balco
                +a1_Be*exp(e1_Be*rho*100*burial(kk))...
                +a2_Be*exp(e2_Be*rho*100*burial(kk))...
                +a3_Be*exp(e3_Be*rho*100*burial(kk)));

          Pz_Al = pfac(kk)*(Pspal_Al*exp(-rho*100*burial(kk)/Lspal_Al)...           % Balco
                +a1_Al*exp(e1_Al*rho*100*burial(kk))...
                +a2_Al*exp(e2_Al*rho*100*burial(kk))...
                +a3_Al*exp(e3_Al*rho*100*burial(kk)));

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
    plot([min(NBe_fin) max(NBe_fin)],[min(NBe_fin) max(NBe_fin)]*R_SLHL,'-r')
    xlabel('^1^0Be [at/g]')
    ylabel('^2^6Al [at/g]')

cd ..\..; close(h); toc

end


function nUpdateWaitbar(~)
    global p h
        t = toc;
        waitbar(p/1000, h, ['Running simulations, please wait... ',num2str(round(t/(p/1000)-t)),' seconds']);
        p = p + 1;
end

