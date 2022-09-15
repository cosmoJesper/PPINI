
clear all; close all; clc;
rng(1)
site = '';

T = readtable('Mohlin.xlsx','VariableNamingRule','preserve');
names = table2array(T(:,1));
NBe = table2array(T(:,2));
NBe_unc = table2array(T(:,3));
NAl = table2array(T(:,4));
NAl_unc = table2array(T(:,5));
depths = table2array(T(:,6));

% NBe(4) = 0;
% NAl(4) = 0;
% NBe_unc(4) = mean(NBe_unc)/5;
% NAl_unc(4) = mean(NAl_unc)/5;
% names{4} = '0';

k = (1:numel(NAl))';
% k = (16:25)';
k = [1 2 3 4 6]';
% k = k';
% d = 450;
% k = find(depths == d);
% k(5) = [];

% lat = 63.86189;                                             % Latitude of sink site
% lon = -148.85488;                                           % Longitude of sink site
% elev = 563;                                        % Elevation of sampling site [m asl.]
% R = 6.87;

lat = 47;                                             % Latitude of sink site
lon = 8;                                           % Longitude of sink site
elev = 379;                                        % Elevation of sampling site [m asl.]
R = 7.6;

TBe = 1.387e6;                                          % 10Be half-life [yr]
TAl = 0.705e6;                                          % 26Al half-life [yr]
lambda_Be = log(2)/TBe;                                 % 10Be mean lifetime [yr^-1]
lambda_Al = log(2)/TAl;                                 % 26Al mean lifetime [yr^-1]


% [age_final,error_final,m,b] = iso_bur_implement(NAl(k),NAl_unc(k),NBe(k),NBe_unc(k),P_Be*6.8,P_Be);
[age_final,error_final,m,m_std,b,b_std] = isochron_mc_LSD_call(NAl(k),NAl_unc(k),NBe(k),NBe_unc(k),lon,lat,elev,R);
ages(1) = age_final*1e6;
errors(1) = error_final*1e6;

% m = 2;
% b = -1000;
% m_std = 0.3;
% b_std = 100;
% NBe = NBe*1e-5; NAl = NAl*1e-6; NBe_unc = NBe_unc*1e-6; NAl_unc = NAl_unc*1e-6; b = b*1e-6; b_std = b_std*1e-6;

% figure('DefaultAxesFontSize',18)
figure(1)
xlimit = [0 max(NBe(k)+NBe_unc(k))*1.05];
xlim(xlimit)
lh(1) = plot(xlimit,R*xlimit,'-k','LineWidth',1);
hold on
lh(2) = plot(xlimit,m*xlimit+b,'-r','LineWidth',1);
y_plus = (m+m_std)*xlimit+(b+b_std);
y_minus = (m-m_std)*xlimit+(b-b_std);
% plot(xlimit,y_plus,'--r','LineWidth',1);
% plot(xlimit,y_minus,'--r','LineWidth',1);
patch([xlimit fliplr(xlimit)],[y_plus fliplr(y_minus)],'r','FaceAlpha',0.1,'LineStyle','none')

t = 0:0.1:2*pi;
for i=k'
    patch(NBe(i)+NBe_unc(i)*cos(t),NAl(i)+NAl_unc(i)*sin(t),[0.9 0.9 0.9],'LineStyle','-','FaceAlpha',0.5,'LineWidth',1.1)
end

% text(0.05,0.95,...
%     {['Isochron age: ',num2str(round(ages(1)*1e-3)),' +/- ',num2str(round(errors(1)*1e-3)),' kyr']},...
%     'FontSize',14,'FontWeight','normal','Units','normalized');

% textfit(NBe(k),NAl(k),names(k),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',8,'Color','k')
% plot(NBe(k),NAl(k),'+k')
plot(NBe,NAl,'+k')
i = 5;
plot(NBe(i)+NBe_unc(i)*cos(t),NAl(i)+NAl_unc(i)*sin(t),'--k','LineWidth',1.1)

title(site)
xlabel('^1^0Be [10^4 at/g]')
ylabel('^2^6Al [10^4 at/g]')

xticks([0 0.5 1 1.5 2 2.5]*1e4)
xticklabels({'0','0.5','1','1.5','2','2.5'})
yticks([0 5 10 15 20]*1e4)
yticklabels({'0','5','10','15','20'})
ax =gca; ax.XLim = [0 2.8e4]; ax.YLim = [0 21e4];

% legend(lh,['P_2_6/P_1_0 = ',num2str(R)],'Isochron fit','location','northwest')

figure(1)
set(gcf,'Units','normalized','position',[0.6 0.5 0.31 0.42]);
box on
% ax =gca; ax.YLim = [0 15e3];
% m=exp(-2.34e6*(lambda_Al-lambda_Be))*6.87
% plot(xlimit,xlimit*m-1500,'-b')
% 
% [a, b, sigma_a, sigma_b, S] = york_fit(NBe(k)',NAl(k)',NBe_unc(k)',NAl_unc(k)')
% plot(xlimit,xlimit*b+a,'-g')

% legend('boxoff')

set(gcf,'Units','normalized','position',[0.05,0.05,0.3,0.4]);


simple_figure()
 h = text(2e4,16.1e4,['^2^6P/^1^0P = ',num2str(R)]);
set(h,'Rotation',39.5);
