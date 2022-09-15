function MC_AlBe_sink_ero

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

close all;
delete(gcp('nocreate'))
cd code

% Name the output source file and determine the Number of simulations you want:
load mc_source_Mohlin.mat                                % Load file with sink simulations
name = 'mc_sink_Mohlin';                                 % Filename for saving

cd ..\samples
T = readtable('Mohlin.xlsx');                            % Excel file with sample data
cd ..\code

lat = 47;                                                      % Latitude of sink site
lon = 8;                                                      % Longitude of sink site
elev = 379;                                           % Elevation of sampling site [m asl.]

%*********** Sink history ************
Nt = 51;                                              % Number of sink steps.
sink_time = linspace(0,0.4,Nt)*1e6;                   % Possible burial times in sink [Myr]
Ne = 51;                                              % Number of erosion steps to make sure all erosionrates are used.
erates = 1*linspace(0,100,Ne)*1e-6;                   % Possible erosion rates in sink [m/yr]

%*********** Other info **************
rho = 2.0;                                            % Density of sediment in the sink [g/cm3]
P_SLHL = 4;
R = 6.75;
dt = 2000;                                            % dt [yr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%           Code form here           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depths = table2array(T(:,6));
depths = 3000;

if -max(sink_time)*min(erates) >= max(depths)/100
    max_accu = max(depths)/100/max(sink_time);
    msg = ['Maximum accumulation rate allowed is ',num2str(max_accu),' m/Myr'];
    error(msg)
end

LSD_DG(lat, lon, elev,0,2e6,-1,10);
load LSDout;
P_time = P_SLHL*mean(LSDout.Be);
Pspal_Be = P_time(1);
LSD_DG(lat, lon, elev,0,2e6,-1,26);
load LSDout;
P_time = R*P_SLHL*mean(LSDout.Al);
Pspal_Al = P_time(1);

% Loading in from previous source run
NBe_fin = mc_source.NBe_fin;
NAl_fin = mc_source.NAl_fin;
Nsim = mc_source.Nsim;
                               
% Creating an unifrom distribution of erosion rates
c = floor(Nsim/Nt);
erate_p = min(erates)*ones(1,c); 
for i = 1:Ne
    erate_p(floor(c/Ne)*(i-1)+1:floor(c/Ne)*i) = erates(i);
end
erate_p(floor(c/Ne)*Ne+1:end) = erates(randi(Ne,1,c-floor(c/Ne)*Ne));

erate_pool = min(erates)*ones(1,Nsim);
for i = 1:Nt
    erate_pool(c*(i-1)+1:c*i) = erate_p;
end
erate_pool(c*Nt+1:end) = erates(randi(Ne,1,Nsim-c*Nt));

%************* Time settings ***************
time = 0:dt:max(sink_time);

%************* Cosmo info **************                
Pnmc_Be = Pspal_Be*0.015;                               % 10Be production by neutron-capture [atoms/(g*yr)]
Pfm_Be = Pspal_Be*0.005;                                % 10Be production by fast muons [atoms/(g*yr)]
P_Be = Pspal_Be + Pnmc_Be + Pfm_Be;                     % Total 10Be production
Pnmc_Al = Pspal_Al*0.018;                               % 26Al production by neutron-capture [atoms/(g*yr)]
Pfm_Al = Pspal_Al*0.006;                                % 26Al production by fast muons [atoms/(g*yr)]
P_Al = Pspal_Al + Pnmc_Al + Pfm_Al;                     % Total 26Al production


gmr = -0.03417; % Assorted constants
dtdz = 0.0065; % Lapse rate from standard atmosphere
pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev.*dtdz)) ) );
load consts_LSD
z = linspace(0,round(max(erates)*max(sink_time)+max(depths)/100+1,-1),1000);

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

P_mu_Al = P_mu_total_edited(z*100*rho,pressure,c26,'no');
P_mu_Be = P_mu_total_edited(z*100*rho,pressure,c10,'no');

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
x0_Be = [1 -rho*100/Lnmc_Be 0.1 -rho*100/Lfm_Be 0.001 -rho*100/Lfm_Be]; 
x0_Al = [1 -rho*100/Lnmc_Al 0.1 -rho*100/Lfm_Al 0.001 -rho*100/Lfm_Al]; 

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
erate_fin = zeros(1,Nsim);   

figure;
    hold on
    Pz_Be1 = Pspal_Be*exp(-rho*100*z/Lspal_Be)...
        +a1_Be*exp(e1_Be*z)...
        +a2_Be*exp(e2_Be*z)...
        +a3_Be*exp(e3_Be*z);

    Pz_Al1 = Pspal_Al*exp(-rho*100*z/Lspal_Al)...
        +a1_Al*exp(e1_Al*z)...
        +a2_Al*exp(e2_Al*z)...
        +a3_Al*exp(e3_Al*z);

    Pz_Be2 = Pspal_Be*exp(-rho*100*z/Lspal_Be)...
        +Pnmc_Be*exp(-rho*100*z/Lnmc_Be)...
        +Pfm_Be*exp(-rho*100*z/Lfm_Be);

    Pz_Al2 = Pspal_Al*exp(-rho*100*z/Lspal_Al)...
        +Pnmc_Al*exp(-rho*100*z/Lnmc_Al)...
        +Pfm_Al*exp(-rho*100*z/Lfm_Al);
    
    subplot(1,2,1)
        hold on
        plot(Pz_Be1,z,'-r')
        plot(P_mu_Be + Pspal_Be*exp(-rho*100*z/Lspal_Be),z,'--k')
        plot(Pz_Be2,z,'-b')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 10Be production [at/g/yr]')
        ylabel('Z [m]')
        legend('Balco2017','Balco2017 fit','Heisinger','location','southeast')
    subplot(1,2,2)
        hold on
        plot(Pz_Al1,z,'-r')
        plot(P_mu_Al + Pspal_Al*exp(-rho*100*z/Lspal_Al),z,'--k')
        plot(Pz_Al2,z,'-b')
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [m]')
        legend('Balco2017','Balco2017 fit','Heisinger','location','southeast')

global h p
parpool;

for d = 1e-2*unique(depths)'
    
    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    afterEach(D, @nUpdateWaitbar);
    p = 1;

    tic

    parfor (j=1:Nsim,3)
        Age = sink_time(randperm(length(sink_time),1));
        time_steps = find(Age >= time);                         % timesteps to calculate
        ntime = length(time_steps);                             % number of timesteps

        erate = erate_pool(j);                                  % picking an erosion rate

        burial = zeros(1,ntime);
        burial(1) = d;

        if Age ~= 0
            for i=2:ntime-1                                 % Looping to calculate burial depth with time (stepping backwards in time)
                    burial(i) = burial(i-1) + dt*erate;      % else normal procedure for interglacial
            end
            burial(ntime) = burial(ntime-1) + dt*erate;
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
            Pz_Be = Pspal_Be*exp(-rho*100*burial(ik)/Lspal_Be)...
                +a1_Be*exp(e1_Be*burial(ik))...
                +a2_Be*exp(e2_Be*burial(ik))...
                +a3_Be*exp(e3_Be*burial(ik));

            Pz_Al = Pspal_Al*exp(-rho*100*burial(ik)/Lspal_Al)...
                +a1_Al*exp(e1_Al*burial(ik))...
                +a2_Al*exp(e2_Al*burial(ik))...
                +a3_Al*exp(e3_Al*burial(ik));
                
%                 Pz_Be = Pspal_Be*exp(-rho*100*burial(ik)/Lspal_Be)...
%                     +Pnmc_Be*exp(-rho*100*burial(ik)/Lnmc_Be)...
%                     +Pfm_Be*exp(-rho*100*burial(ik)/Lfm_Be);
% 
%                 Pz_Al = Pspal_Al*exp(-rho*100*burial(ik)/Lspal_Al)...
%                     +Pnmc_Al*exp(-rho*100*burial(ik)/Lnmc_Al)...
%                     +Pfm_Al*exp(-rho*100*burial(ik)/Lfm_Al);
                
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
        erate_fin(j) = erate;  
        
        
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
    
    rate = erate_fin;

    save([name,'_test'],'NBe_inh_initial','NAl_inh_initial','NBe_post',...
        'NAl_post','NBe_tot','NAl_tot','rate','Ages_fin','storage_time','storage_depth',...
        'source_pluck','Tice','erate_source','Nt','Ne','elevation_source','-v7.3')

    %%%%%%%%%%%%%%%% Creating figures %%%%%%%%%%%%%%%%%%%%%
    figure
        title('1:1000 of simulations')
        plot(NBe_tot(1:1000:end),NAl_tot(1:1000:end),'.k')
        hold on
        plot([min(NBe_tot) max(NBe_tot)],[min(NBe_tot) max(NBe_tot)]*7,'-r')

    close(h)

end
delete(gcp('nocreate'));

end

function nUpdateWaitbar(~)
    global p h
        t = toc;
        waitbar(p/1000, h, ['Estimated waiting time: ',num2str(round(t/(p/1000)-t)),' seconds']);
        p = p + 1;
end
