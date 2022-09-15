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
cd other/code                                  % Change directory
load mc_source_Mohlin.mat                        % Load file with sink simulations
name = 'mc_sink_Mohlin';                                 % Filename for saving
file = 'Mohlin.xlsx';

%*********** Geographical settings ************

lat = 47;                                              % Latitude of sink site
lon = 8;                                            % Longitude of sink site
elev = 379;                                              % Elevation of sampling site [m asl.]

%*********** Sink history ************
min_sink_time = 0;                                    % Possible burial times in sink [Myr]         
max_sink_time = 0.3;                                   % Possible burial times in sink [Myr]
min_accum_rate = 0;                                   % Possible erosion rates in sink [g/cm2/Myr]
max_accum_rate = 6000;                               % Possible erosion rates in sink [g/cm2/Myr]
starting_depth = 6000;                                % Estimated starting depth of sample [g/cm2]

%*********** Other info **************
P_SLHL = 4;                                           % 10Be production sea-level high-latitude [at/g/yr]
R = 6.75;                                             % 26Al/10Be production ratio at sea-level high-latitude

dt = 2000;                                            % dt [yr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           Code form here           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
delete(gcp('nocreate'))
cd ..\samples
T = readtable(file);                            % Excel file with sample data
cd ..\code

Nt = 51;                                              % Number of sink steps.
sink_time = linspace(min_sink_time,max_sink_time,Nt)*1e6;                    % Possible burial times in sink [Myr]
Ne = Nt;                                              % Number of erosion steps to make sure all erosionrates are used.
accums = linspace(min_accum_rate,max_accum_rate,Ne)*1e-6;                 % Possible erosion rates in sink [g/cm2/yr]

LSD_DG(lat, lon, elev,0,2e6,-1,10);
load LSDout;
P_time = P_SLHL*mean(LSDout.Be);
Pspal_Be = P_time(1);                       % 10Be production by spallation [atoms/(g*yr)]
LSD_DG(lat, lon, elev,0,2e6,-1,26);
load LSDout;
P_time = R*P_SLHL*mean(LSDout.Al);
Pspal_Al = P_time(1);                       % 26Al production by spallation [atoms/(g*yr)]

% Loading in from previous source run
NBe_fin = mc_source.NBe_fin;
NAl_fin = mc_source.NAl_fin;
Nsim = mc_source.Nsim;

%************* Time settings ***************
time = 0:dt:max(sink_time);

%************* Cosmo info **************                
                                   
Pnmc_Be = Pspal_Be*0.015;                               % 10Be production by neutron-capture [atoms/(g*yr)]
Pfm_Be = Pspal_Be*0.005;                                % 10Be production by fast muons [atoms/(g*yr)]
Pnmc_Al = Pspal_Al*0.018;                               % 26Al production by neutron-capture [atoms/(g*yr)]
Pfm_Al = Pspal_Al*0.006;                                % 26Al production by fast muons [atoms/(g*yr)]

gmr = -0.03417; % Assorted constants
dtdz = 0.0065; % Lapse rate from standard atmosphere
pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev.*dtdz)) ) );
load consts_LSD
z = linspace(0,round(max(accums)*max(sink_time)+starting_depth,-1),1000);

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

P_mu_Al = P_mu_total_edited(z,pressure,c26,'no');
P_mu_Be = P_mu_total_edited(z,pressure,c10,'no');

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

fitfun = fittype( @(a1,e1,a2,e2,a3,e3,x) a1*exp(x*e1) + a2*exp(x*e2) + a3*exp(x*e3));
x0_Be = [1 -1/Lnmc_Be 0.1 -1/Lfm_Be 0.001 -1/Lfm_Be]; 
x0_Al = [1 -1/Lnmc_Al 0.1 -1/Lfm_Al 0.001 -1/Lfm_Al]; 

options = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0 -10 0 -10 0 -10], 'Upper',[10 0 10 0 10 0],'StartPoint',x0_Be);
                 
[fitted_curve_Be,gof] = fit(z',P_mu_Be',fitfun,options);
a1_Be = fitted_curve_Be.a1;
e1_Be = fitted_curve_Be.e1;
a2_Be = fitted_curve_Be.a2;
e2_Be = fitted_curve_Be.e2;
a3_Be = fitted_curve_Be.a3;
e3_Be = fitted_curve_Be.e3;

options = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0 -10 0 -10 0 -10], 'Upper',[10 0 10 0 10 0],'StartPoint',x0_Al);
                 
[fitted_curve_Al,gof] = fit(z',P_mu_Al',fitfun,options);
a1_Al = fitted_curve_Al.a1;
e1_Al = fitted_curve_Al.e1;
a2_Al = fitted_curve_Al.a2;
e2_Al = fitted_curve_Al.e2;
a3_Al = fitted_curve_Al.a3;
e3_Al = fitted_curve_Al.e3;

NBe_inh_initial = zeros(1,Nsim);                        % Lists for saving
NAl_inh_initial = zeros(1,Nsim);    
NBe_post = zeros(1,Nsim);
NAl_post = zeros(1,Nsim);
NBe_tot = zeros(1,Nsim);
NAl_tot = zeros(1,Nsim); 
Ages_fin = zeros(1,Nsim); 
accum_fin = zeros(1,Nsim);   

figure;
    hold on
    Pz_Be1 = Pspal_Be*exp(-z/Lspal_Be)...
        +a1_Be*exp(e1_Be*z)...
        +a2_Be*exp(e2_Be*z)...
        +a3_Be*exp(e3_Be*z);

    Pz_Al1 = Pspal_Al*exp(-z/Lspal_Al)...
        +a1_Al*exp(e1_Al*z)...
        +a2_Al*exp(e2_Al*z)...
        +a3_Al*exp(e3_Al*z);

    Pz_Be2 = Pspal_Be*exp(-z/Lspal_Be)...
        +Pnmc_Be*exp(-z/Lnmc_Be)...
        +Pfm_Be*exp(-z/Lfm_Be);

    Pz_Al2 = Pspal_Al*exp(-z/Lspal_Al)...
        +Pnmc_Al*exp(-z/Lnmc_Al)...
        +Pfm_Al*exp(-z/Lfm_Al);
    
    subplot(1,2,1)
        hold on
        plot(Pz_Be1,z,'-r')
        plot(P_mu_Be + Pspal_Be*exp(-z/Lspal_Be),z,'--k')
        plot(Pz_Be2,z,'-b')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 10Be production [at/g/yr]')
        ylabel('Z [g/cm2]]')
        legend('Balco2017','Balco2017 fit','Heisinger','location','southeast')
    subplot(1,2,2)
        hold on
        plot(Pz_Al1,z,'-r')
        plot(P_mu_Al + Pspal_Al*exp(-z/Lspal_Al),z,'--k')
        plot(Pz_Al2,z,'-b')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [g/cm2]')
        legend('Balco2017','Balco2017 fit','Heisinger','location','southeast')


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

        accum = accums(randperm(length(accums),1));             % picking an erosion rate

        burial = zeros(1,ntime);
        burial(1) = starting_depth + accum*Age;

        if Age ~= 0
            for i=2:ntime-1                                 % Looping to calculate burial depth with time (stepping backwards in time)
                    burial(i) = burial(i-1) - dt*accum;      % else normal procedure for interglacial
            end
            burial(ntime) = burial(ntime-1) - dt*accum;
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
            Pz_Be = Pspal_Be*exp(-burial(ik)/Lspal_Be)...           % Balco
                +a1_Be*exp(e1_Be*burial(ik))...
                +a2_Be*exp(e2_Be*burial(ik))...
                +a3_Be*exp(e3_Be*burial(ik));

            Pz_Al = Pspal_Al*exp(-burial(ik)/Lspal_Al)...           % Balco
                +a1_Al*exp(e1_Al*burial(ik))...
                +a2_Al*exp(e2_Al*burial(ik))...
                +a3_Al*exp(e3_Al*burial(ik)); 
           
%                 Pz_Be = Pspal_Be*exp(-burial(ik)/Lspal_Be)...     % Heisinger
%                     +Pnmc_Be*exp(-burial(ik)/Lnmc_Be)...
%                     +Pfm_Be*exp(-burial(ik)/Lfm_Be);
% 
%                 Pz_Al = Pspal_Al*exp(-burial(ik)/Lspal_Al)...   % Heisinger
%                     +Pnmc_Al*exp(-burial(ik)/Lnmc_Al)...
%                     +Pfm_Al*exp(-burial(ik)/Lfm_Al);
                
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
        accum_fin(j) = accum;  
        
        
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

    rate = accum_fin;

    save([name,'_test'],'NBe_inh_initial','NAl_inh_initial','NBe_post',...
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

    close(h)

delete(gcp('nocreate'));
cd ..\..
end

function nUpdateWaitbar(~)
    global p h
        t = toc;
        waitbar(p/1000, h, ['Estimated waiting time: ',num2str(round(t/(p/1000)-t)),' seconds']);
        p = p + 1;
end
