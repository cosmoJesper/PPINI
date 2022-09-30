function MC_AlBe_sink_accum

% Sink sample generator for the P-PINI code.
%
% The code will generate a .mat file containing:
%   10Be concentrations [atoms/gr]
%   26Al concentrations [atoms/gr]
%   List of input parameters
%
% The workflow of the program is:
%
%         Run MC_AlBe_source.m
%                   |
% ->       Run MC_AlBe_sink.m
%                   |
%         Run PPINI_analylsis.m

% Name the output source file and determine the Number of simulations you want:
cd other/code                                           % Change directory
load mc_source_Mohlin.mat                               % Load file with sink simulations
name = 'mc_sink_Mohlin';                                % Filename for saving
file = 'Mohlin.xlsx';                                  % Name of excelfile in sample folder

%********** Accumulation or erosion ************
ero = 0;  % Enter 1 for erosion, 0 for accumulation

%*********** Geographical settings ************
lat = 47;                                             % Latitude of sink site
lon = 8;                                           % Longitude of sink site
elev = 379;                                             % Elevation of sampling site [m asl.]

%*********** Sink history ************
min_sink_time = 0;                                      % Possible burial times in sink [Myr]         
max_sink_time = 0.3;                                      % Possible burial times in sink [Myr]
min_rate = 0;                                           % Possible erosion/accumulation rates in sink [g/cm2/Myr] (Depending on ero setting)
max_rate = 6000;                                       % Possible erosion/accumulation rates in sink [g/cm2/Myr]

depth = 6000;                                         % IF ACCUMULATION: Input estimated depth of samples within unit [g/cm2]
                                                        % IF EROSION: Input estimated sampling depth of samples [g/cm2]
%*********** Other info **************
P_SLHL_Be = 4;                                        % 10Be production sea-level high-latitude [at/g/yr]
R_SLHL = 6.8;                                         % 26Al/10Be production ratio at SLHL
P_SLHL_Al = R_SLHL*P_SLHL_Be;                         % 26Al production at SLHL [at/g/yr]

dt = 2000;                                            % dt [yr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           Code form here           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; delete(gcp('nocreate')); cd ..\samples
T = readtable(file);                                                    % Excel file with sample data
cd ..\code

Nt = 51;                                                                % Number of sink steps.
sink_time = linspace(min_sink_time,max_sink_time,Nt)*1e6;               % Possible burial times in sink [Myr]
Ne = Nt;                                                                % Number of erosion steps
rates = linspace(min_rate,max_rate,Ne)*1e-6;               % Possible erosion rates in sink [g/cm2/yr]

LSD_DG(lat, lon, elev,0,2e6,-1,10);
load LSDout;
P_time = P_SLHL_Be*mean(LSDout.Be);
Pspal_Be = P_time(1);                       % 10Be production by spallation [atoms/(g*yr)]
LSD_DG(lat, lon, elev,0,2e6,-1,26);
load LSDout;
P_time = P_SLHL_Al*mean(LSDout.Al);
Pspal_Al = P_time(1);                       % 26Al production by spallation [atoms/(g*yr)]

% Loading in from previous source run
NBe_fin = mc_source.NBe_fin;
NAl_fin = mc_source.NAl_fin;
Nsim = mc_source.Nsim;

%************* Setting up ***************
time = 0:dt:max(sink_time);
z = linspace(0,min([1.9e5 round(max(rates)*max(sink_time)+depth,-1)]),1000);

%************* Cosmo info **************                
gmr = -0.03417;                                         % Assorted constants
dtdz = 0.0065;                                          % Lapse rate from standard atmosphere
pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev.*dtdz)) ) ); % Calculating atmospheric pressure at site

c10.fC = 0.704;
c10.fD = 0.1828;
c10.Natoms = 2.006e22;
c10.fstar = 0.00157;
c10.sigma190 = 37.8e-30;
c26.fC = 0.296;
c26.fD = 0.6559;
c26.Natoms = 1.003e22;
c26.fstar = 0.0118;
c26.sigma190 = 521e-30;

P_mu_Be = P_mu_total_edited(z,pressure,c10,'no');
P_mu_Al = P_mu_total_edited(z,pressure,c26,'no');

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

objfcn = @(b,x) b(1)*exp(x*b(4)) + b(2)*exp(x*b(5)) + b(3)*exp(x*b(6));
nlm_Be = fitnlm(z',P_mu_Be',objfcn,[1 1 1 -1/Lnmc_Be -1/Lfm_Be -1/Lfm_Be]);
b_Be = nlm_Be.Coefficients.Estimate;
nlm_Al = fitnlm(z',P_mu_Al',objfcn,[1 1 1 -1/Lnmc_Al -1/Lfm_Al -1/Lfm_Al]);
b_Al = nlm_Al.Coefficients.Estimate;

NBe_inh_initial = zeros(1,Nsim);                        % Lists for saving
NAl_inh_initial = zeros(1,Nsim);    
NBe_post = zeros(1,Nsim);
NAl_post = zeros(1,Nsim);
NBe_tot = zeros(1,Nsim);
NAl_tot = zeros(1,Nsim); 
Ages_fin = zeros(1,Nsim); 
rate_fin = zeros(1,Nsim);   

figure;
    hold on
    Pz_Be1 = Pspal_Be*exp(-z/Lspal_Be)...
        +b_Be(1)*exp(b_Be(4)*z)...
        +b_Be(2)*exp(b_Be(5)*z)...
        +b_Be(3)*exp(b_Be(6)*z);

    Pz_Al1 = Pspal_Al*exp(-z/Lspal_Al)...
        +b_Al(1)*exp(b_Al(4)*z)...
        +b_Al(2)*exp(b_Al(5)*z)...
        +b_Al(3)*exp(b_Al(6)*z);

    subplot(1,3,1)
        hold on
        plot(P_mu_Be + Pspal_Be*exp(-z/Lspal_Be),z,'-r')
        plot(Pz_Be1,z,'--k')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 10Be production [at/g/yr]')
        ylabel('Z [g/cm2]]')
        legend('Balco2017','Balco2017 fit','location','southeast')
    subplot(1,3,2)
        hold on
        plot(P_mu_Al + Pspal_Al*exp(-z/Lspal_Al),z,'-r')
        plot(Pz_Al1,z,'--k')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [g/cm2]')
        legend('Balco2017','Balco2017 fit','location','southeast')
    subplot(1,3,3)
        hold on
        plot((P_mu_Al + Pspal_Al*exp(-z/Lspal_Al))./(P_mu_Be + Pspal_Be*exp(-z/Lspal_Be)),z,'-r')
        plot(Pz_Al1./Pz_Be1,z,'--k')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [g/cm2]')
        legend('Balco2017','Balco2017 fit','location','southeast')

global h p
parpool;

    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    afterEach(D, @nUpdateWaitbar);
    p = 1;

    tic

    parfor (j=1:Nsim,3)
        Age = sink_time(randperm(length(sink_time),1));
        time_steps = find(Age >= time);                         % timesteps to calculate
        ntime = length(time_steps);                             % number of timesteps

        rate = rates(randperm(length(rates),1));             % picking an erosion/accumulation rate

        burial = zeros(1,ntime);
        if ero == 1
            burial(1) = depth;
            if Age ~= 0
                for i=2:ntime-1                                 % Looping to calculate burial depth with time (stepping backwards in time)
                        burial(i) = burial(i-1) + dt*rate;      % else normal procedure for interglacial
                end
                burial(ntime) = burial(ntime-1) + dt*rate;
            end
        else
            burial(1) = depth + rate*Age;
            if Age ~= 0
                for i=2:ntime-1                                 % Looping to calculate burial depth with time (stepping backwards in time)
                        burial(i) = burial(i-1) - dt*rate;      % else normal procedure for interglacial
                end
                burial(ntime) = burial(ntime-1) - dt*rate;
            end
        end
        
        NBe_inh_run = zeros(ntime,1);
        NAl_inh_run = zeros(ntime,1);
        NBe_post_run = zeros(ntime,1);
        NAl_post_run = zeros(ntime,1);

        NBe_inh_run(ntime) = NBe_fin(j);
        NAl_inh_run(ntime) = NAl_fin(j);   
        NBe_post_run(ntime) = 0;
        NAl_post_run(ntime) = 0;

        for ik = ntime-1:-1:1

            Pz_Be = b_Be(1)*exp(burial(ik)*b_Be(4))...
                + b_Be(2)*exp(burial(ik)*b_Be(5))...
                + b_Be(3)*exp(burial(ik)*b_Be(6))...
                + Pspal_Be*exp(-burial(ik)/Lspal_Be);
            Pz_Al = b_Al(1)*exp(burial(ik)*b_Al(4))...
                + b_Al(2)*exp(burial(ik)*b_Al(5))...
                + b_Al(3)*exp(burial(ik)*b_Al(6))...
                + Pspal_Al*exp(-burial(ik)/Lspal_Al);

            NBe_inh_run(ik) = NBe_inh_run(ik+1)*exp(-lambda_Be*dt);
            NAl_inh_run(ik) = NAl_inh_run(ik+1)*exp(-lambda_Al*dt);
            NBe_post_run(ik) = NBe_post_run(ik+1)*exp(-lambda_Be*dt) + Pz_Be*(1-exp(-dt*lambda_Be))/lambda_Be;
            NAl_post_run(ik) = NAl_post_run(ik+1)*exp(-lambda_Al*dt) + Pz_Al*(1-exp(-dt*lambda_Al))/lambda_Al;        
        end

        NBe_inh_initial(j) = NBe_fin(j);
        NAl_inh_initial(j) = NAl_fin(j);    
        NBe_post(j) = NBe_post_run(1);
        NAl_post(j) = NAl_post_run(1);
        NBe_tot(j) = NBe_inh_run(1) + NBe_post_run(1);
        NAl_tot(j) = NAl_inh_run(1) + NAl_post_run(1); 
        Ages_fin(j) = Age;
        rate_fin(j) = rate;  

        if j==(Nsim*0.001)*floor(j/(Nsim*0.001)) 
            send(D, j);                             % Updating waitbar
        end    
    end

    waitbar(1,h,'Saving...')

    %%%%%%%%%%%%%%%%% Saving results %%%%%%%%%%%%%%%%%%%%%%

    storage_time = mc_source.storage_time;
    storage_depth = mc_source.storage_depth;
    Tice = mc_source.Tice;  
    erate_source = mc_source.erate_source;
    elevation_source = mc_source.elevation_source;
    source_pluck = mc_source.source_pluck;

    rate = rate_fin;

    save(name,'NBe_inh_initial','NAl_inh_initial','NBe_post',...
        'NAl_post','NBe_tot','NAl_tot','rate','Ages_fin','storage_time','storage_depth',...
        'source_pluck','Tice','erate_source','Nt','Ne','elevation_source','-v7.3')

    %%%%%%%%%%%%%%%% Creating figures %%%%%%%%%%%%%%%%%%%%%
    figure
        plot(NBe_tot(1:1000:end),NAl_tot(1:1000:end),'.k')
        hold on
        title('1:1000 of sink simulations')
        plot([min(NBe_tot) max(NBe_tot)],[min(NBe_tot) max(NBe_tot)]*7,'-r')
        xlabel('^1^0Be [at/g]')
        ylabel('^2^6Al [at/g]')

    close(h); delete(gcp('nocreate')); cd ..\..

end

function nUpdateWaitbar(~)
    global p h
        t = toc;
        waitbar(p/1000, h, ['Estimated waiting time: ',num2str(round(t/(p/1000)-t)),' seconds']);
        p = p + 1;
end
