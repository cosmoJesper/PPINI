close all; clear all; clc;

lat = 36;                                                      % Latitude of sink site
lon = 81.5;                                                      % Longitude of sink site
elev = 5708;                                         % Elevation of sampling site [m asl.]
rho = 2.6;                                            % Density of sediment in the sink [g/cm3]
P_SLHL = 4;
R = 6.8;

LSD_DG(lat, lon, elev,0,2e7,-1,10);
load LSDout;
mbe = mean(LSDout.Be);
P_time = P_SLHL*mean(LSDout.Be);
Pspal_Be = P_time(1);
LSD_DG(lat, lon, elev,0,2e7,-1,26);
load LSDout;
mal = mean(LSDout.Al);
P_time = R*P_SLHL*mean(LSDout.Al);
Pspal_Al = P_time(1);

                                   % 10Be production by spallation [atoms/(g*yr)]
Pnmc_Be = Pspal_Be*0.015;                               % 10Be production by neutron-capture [atoms/(g*yr)]
Pfm_Be = Pspal_Be*0.005;                                % 10Be production by fast muons [atoms/(g*yr)]
P_Be = Pspal_Be + Pnmc_Be + Pfm_Be;                     % Total 10Be production
% Pspal_Al = Pspal_Be*R;                                % 26Al production by spallation [atoms/(g*yr)]
Pnmc_Al = Pspal_Al*0.018;                               % 26Al production by neutron-capture [atoms/(g*yr)]
Pfm_Al = Pspal_Al*0.006;                                % 26Al production by fast muons [atoms/(g*yr)]
P_Al = Pspal_Al + Pnmc_Al + Pfm_Al;      

gmr = -0.03417; % Assorted constants
dtdz = 0.0065; % Lapse rate from standard atmosphere
% pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (elev.*dtdz)) ) );
pressure = 989.1 .* exp( elev./(-7588));
% pressure = ERA40atm(lat,lon,elev);
load consts_LSD
Sphi = consts.SPhiInf;
Rc = 0;

f10C = 0.704;
f10D = 0.1828;
c10.fC = f10C;
c10.fD = f10D;

f26C = 0.296;
f26D = 0.6559;
c26.fC = f26C;
c26.fD = f26D;

z = linspace(0,300,1000);
c26.Natoms = consts.Natoms26;
% c26.k_neg = consts.k_neg26;
c26.k_neg = 0.0118;
% c26.sigma190 = consts.sigma190_26;
c26.sigma190 = 521e-30;
c26.mfluxRef = consts.mfluxRef;
P_mu_Al = P_mu_total_full(z*100*rho,pressure,c26,'no');
% P_mu_Al = P_mu_totalLSD(z*100*rho,pressure,Rc,Sphi,c26,'no');
c10.Natoms = consts.Natoms10;
% c10.k_neg = consts.k_neg10;
c10.k_neg = 0.00157;
% c10.sigma190 = consts.sigma190_10;
c10.sigma190 = 37.8e-30;
c10.mfluxRef = consts.mfluxRef;
P_mu_Be = P_mu_total_full(z*100*rho,pressure,c10,'no');
% P_mu_Be = P_mu_totalLSD(z*100*rho,pressure,Rc,Sphi,c10,'no');

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

% x0 = [1 1 -rho*100/Lnmc_Be -rho*100/Lfm_Be]; 
% fitfun = fittype( @(a1,a2,e1,e2,x) a1*exp(x*e1) + a2*exp(x*e2));
% [fitted_curve,gof] = fit(z',P_mu_Be',fitfun,'StartPoint',x0);
% P = Pspal_Be*exp(-rho*100*z'/Lspal_Be) + fitted_curve(z');


% x0 = [1 -rho*100/Lnmc_Be 1 -rho*100/Lfm_Be]; 
% % fitfun = fittype( @(a1,a2,e1,e2,x) a1*exp(x*e1) + a2*exp(x*e2) + a3*exp(x*e3));
% [fitted_curve,gof] = fit(z',P_mu_Be','exp2','StartPoint',x0);
% P_fit_Be = Pspal_Be*exp(-rho*100*z'/Lspal_Be) + fitted_curve(z');
% 
% fitted_curve
% 
% x0 = [1 -rho*100/Lnmc_Al 1 -rho*100/Lfm_Al]; 
% % fitfun = fittype( @(a1,a2,e1,e2,x) a1*exp(x*e1) + a2*exp(x*e2) + a3*exp(x*e3));
% [fitted_curve,gof] = fit(z',P_mu_Al','exp2','StartPoint',x0);
% P_fit_Al = Pspal_Al*exp(-rho*100*z'/Lspal_Al) + fitted_curve(z');



% options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
% x0 = [1 -rho*100/Lnmc_Al 0.1 -rho*100/Lfm_Al 0.01 -rho*100/Lfm_Al]; 
% fun = @(x) exp_fmin_func(x,z',P_mu_Al');
% [x,fval] = fminsearch(fun,x0,options)
% P_fit_Al = Pspal_Al*exp(-rho*100*z'/Lspal_Al) + x(1)*exp(x(2)*z') + x(3)*exp(x(4)*z') + x(5)*exp(x(6)*z');

% 
% 
% fitfun = fittype( @(a1,e1,a2,e2,a3,e3,x) a1*exp(x*e1) + a2*exp(x*e2) + a3*exp(x*e3));
% x0 = [1 -rho*100/Lnmc_Be 0.1 -rho*100/Lfm_Be 0.01 -rho*100/Lfm_Be]; 
% [fitted_curve_Be,gof] = fit(z',P_mu_Be',fitfun,'StartPoint',x0);
% P_fit_Be = Pspal_Be*exp(-rho*100*z'/Lspal_Be) + fitted_curve_Be(z');
% x0 = [1 -rho*100/Lnmc_Al 0.1 -rho*100/Lfm_Al 0.01 -rho*100/Lfm_Al]; 
% [fitted_curve_Al,gof] = fit(z',P_mu_Al',fitfun,'StartPoint',x0);
% P_fit_Al = Pspal_Al*exp(-rho*100*z'/Lspal_Al) + fitted_curve_Al(z');

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
P_fit_Be = Pspal_Be*exp(-rho*100*z'/Lspal_Be) + fitted_curve_Be(z');
P_fit_Al = Pspal_Al*exp(-rho*100*z'/Lspal_Al) + fitted_curve_Al(z');

PBe_tot = Pspal_Be*exp(-rho*100*z/Lspal_Be)+P_mu_Be;
PAl_tot = Pspal_Al*exp(-rho*100*z/Lspal_Al)+P_mu_Al;


figure
    subplot(1,3,1)
    hold on
        Pz_Be = Pspal_Be*exp(-rho*100*z/Lspal_Be)...
            +Pnmc_Be*exp(-rho*100*z/Lnmc_Be)...
            +Pfm_Be*exp(-rho*100*z/Lfm_Be);
        plot(Pz_Be,z,'-b')
        hold on
%         Pz_Be_Braucher = Pspal_Be*exp(-rho*100*z/Lspal_Be)...
%             +0.0138*exp(-rho*100*z/Lnmc_Be)...
%             +0.04300*exp(-rho*100*z/Lfm_Be);
        Pz_Be_Braucher = Pspal_Be*exp(-rho*100*z/Lspal_Be)...
           +Pspal_Be*0.0027*exp(-rho*100*z/Lnmc_Be)...
           +Pspal_Be*0.0087*exp(-rho*100*z/Lfm_Be);
        plot(Pz_Be_Braucher,z,'g')
        plot(PBe_tot,z,'-r');
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 10Be production [at/g/yr]')
        ylabel('Z [m]')
        legend('Heisinger2002','Braucher2011','Balco2017 ModelA','location','southeast')
        ax1 = gca;
        
        plot(P_fit_Be,z,'--k')
    subplot(1,3,2)
    hold on
        Pz_Al = Pspal_Al*exp(-rho*100*z/Lspal_Al)...
            +Pnmc_Al*exp(-rho*100*z/Lnmc_Al)...
            +Pfm_Al*exp(-rho*100*z/Lfm_Al);
        plot(Pz_Al,z,'-b')
        hold on
        Pz_Al_Braucher = Pspal_Al*exp(-rho*100*z/Lspal_Al)...
            +Pspal_Al*0.0275*exp(-rho*100*z/Lnmc_Al)...
            +Pspal_Al*2.65e-4*exp(-rho*100*z/Lfm_Al);
        plot(Pz_Al_Braucher,z,'g')
        plot(PAl_tot,z,'-r');
        set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [m]')
        legend('Heisinger2002','Braucher2011','Balco2017','Balco fit','location','southeast')
        ax2 = gca;
        
        plot(P_fit_Al,z,'--k')
    subplot(1,3,3)
        plot(Pz_Al./Pz_Be,z,'-b')
        hold on
        plot(Pz_Al_Braucher./Pz_Be_Braucher,z,'-g')
        plot(PAl_tot./PBe_tot,z,'-r')
        set(gca,'YDir','reverse')
        xlabel('^2^6Al/^1^0Be')
        ylabel('Z [m]')
        legend('Heisinger2002','Braucher2011','Balco2017','Balco fit','location','southeast')
        ax3 = gca;
        
        plot(P_fit_Al./P_fit_Be,z,'--k')
        
        linkaxes([ax1 ax2 ax3],'y')
% x = 60;        
% P_mu_Al(x)./P_mu_Be(x)  
% P_mu_Al(x)
% P_mu_Be(x)

figure
subplot(1,2,1)
    hold on
        plot(P_mu_Be,z,'-r');
        plot(fitted_curve_Be(z'),z,'--b')
%         plot(x(1)*exp(x(2)*z') + x(3)*exp(x(4)*z') + x(5)*exp(x(6)*z'),z,'--k')
%         set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 10Be production [at/g/yr]')
        ylabel('Z [m]')
        legend('Balco2017 ModelA','fit','location','southeast')
subplot(1,2,2)
    hold on
        plot(P_mu_Al,z,'-r');
        plot(fitted_curve_Al(z'),z,'--b')
%         plot(x(1)*exp(x(2)*z') + x(3)*exp(x(4)*z') + x(5)*exp(x(6)*z'),z,'--k')
%         set(gca,'XScale','log')
        set(gca,'YDir','reverse')
        xlabel('Total 26Al production [at/g/yr]')
        ylabel('Z [m]')
        legend('Balco2017 ModelA','fit','location','southeast')
        
      
