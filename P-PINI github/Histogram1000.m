
clear all; close all; clc;
rng(3)

test_age = 1380000; % [yr]
test_erate = 1e-5; % [m/yr]
n_runs = 10;
g_noise = 0.05;
n_samples = 8;

%%% Loading the 'library'
% load mc_sink3_Tromsberg_15m
m = matfile('mc_sink_PuluB2_test');

NBe_post = m.NBe_post;
NAl_post = m.NAl_post;
NBe_tot = m.NBe_tot;
NAl_tot = m.NAl_tot;
erate_fin = m.erate_fin;
Ages_fin = m.Ages_fin;
storage_time = m.storage_time;
storage_depth = m.storage_depth;
erate_source = m.erate_source;
Tice = m.Tice;
source_pluck = m.source_pluck;
elevation_source = m.elevation_source;
NBe_inh = m.NBe_inh_initial;
NAl_inh = m.NAl_inh_initial;
            
R = 6.8;

Pspal_Be = 7.58; %atoms/(g*yr)
Pnmc_Be = Pspal_Be*0.015;
Pfm_Be = Pspal_Be*0.005;
P_Be = Pspal_Be + Pnmc_Be + Pfm_Be;

Pspal_Al = Pspal_Be*R;
Pnmc_Al = Pspal_Al*0.018;
Pfm_Al = Pspal_Al*0.006;
P_Al = Pspal_Al + Pnmc_Al + Pfm_Al;


%%% edges is usefull for plotting histograms
ut = unique(Ages_fin);
dt = ut(3) - ut(2);
edges = linspace(ut(1) - dt/2,ut + dt/2,20);

%%% All the possible burial ages
tburial = ut;

% A = zeros(numel(ub),numel(tburial),numel(k));
ue = unique(erate_fin);

toc_t = 0;

tic

iso_ages = zeros(1,n_runs);
iso_errors = zeros(1,n_runs);
iso_MSWD = zeros(1,n_runs);
mc_ages = zeros(1,n_runs);
mc_errors = zeros(1,n_runs);

[a,idx] = min(abs(ue-test_erate));
erate_samples = ue(idx);
tburial1 = tburial/1e6;
        
% load mc_sink3_Tromsberg_15m
    trim = find(Ages_fin == test_age & NBe_tot < 6e5);
% trim = find(mc_sink.Ages_fin == test_age);
NBe_noise = NBe_tot(trim).*normrnd(1,g_noise,1,length(NBe_tot(trim)));
NAl_noise = NAl_tot(trim).*normrnd(1,g_noise,1,length(NAl_tot(trim)));
erate_noise = mc_sink.erate_fin(trim);
k = 1:n_samples;
filter1 = find(erate_noise == erate_samples);

for n = 1:n_runs

    rand_samp = randperm(numel(filter1),n_samples);
    NAl = NAl_noise(filter1(rand_samp));
    NBe = NBe_noise(filter1(rand_samp));
    NAl_unc = g_noise.*NAl;
    NBe_unc = g_noise.*NBe;
    
    %%% Remove concentrations with too high post-burial production compared to
    %%% the samples.
    [max_conc,Be_idx] = min(NBe);
    max_post_Be = NBe(Be_idx) + NBe_unc(Be_idx); % maximum post burial production we can allow
    [max_conc,Al_idx] = min(NAl);
    max_post_Al = NAl(Al_idx) + NAl_unc(Al_idx);

    to_keep = find(NBe_post < max_post_Be & NAl_post < max_post_Al);

    NBe_tot_run = NBe_tot(to_keep);
    NAl_tot_run = NAl_tot(to_keep);
    erate_fin_run = erate_fin(to_keep);
    Ages_fin_run = Ages_fin(to_keep);

    A = zeros(numel(ue),numel(tburial),numel(k));

    for i = k
        L_NBe{i} = zeros(size(NBe_tot_run));
        L_NAl{i} = zeros(size(NAl_tot_run));
        L{i} = zeros(size(NAl_tot_run));
    end

    for i = k

        % Calculates the likelihood at each simulated data point
        NBe_dd = NBe_tot_run - NBe(i);
        L_NBe{i} = ((2*pi)^2*abs(NBe_unc(i)))^-.5.*exp(-.5*NBe_dd.^2./(NBe_unc(i)).^2);

        NAl_dd= NAl_tot_run - NAl(i);
        L_NAl{i} = ((2*pi)^2*abs(NAl_unc(i)))^-.5.*exp(-.5*NAl_dd.^2./(NAl_unc(i)).^2);

        L{i} = L_NBe{i}.*L_NAl{i};

        maxL = max(L{i});

        % Pick points through rejection sampling
        Pac = L{i}/maxL;
        r = rand(size(L{i}));
        filter = find(r < Pac);

        Ages = Ages_fin_run(filter); 
        erates = erate_fin_run(filter);

        for c1 = 1:numel(ue)
            for c2 = 1:numel(tburial)
                check = find(erates == ue(c1) & Ages == tburial(c2));
                A(c1,c2,i) = numel(check);
            end
        end
    end % current samples

%     [age_final,error_final] = iso_bur_implement(NAl,NAl_unc,NBe,NBe_unc,P_Al,P_Be);
    [age_final,error_final,m,b,MSWD] = isochron_mc_LSD_call(NAl',NAl_unc',NBe',NBe_unc',P_Be,R);
    iso_ages(n) = age_final*1e6;
    iso_errors(n) = error_final*1e6;
    iso_MSWD(n) = MSWD;
    
    A_prod = ones(numel(ue),numel(tburial));

    for i = k
        A_prod = A_prod.*(A(:,:,i)/sum(sum(A(:,:,i))));
    end
    A_prod = A_prod./sum(sum(A_prod));

    A_prod2 = (sum(A_prod)./sum(sum(A_prod)))';
    A_prod2(isnan(A_prod2)) = 0;

    f = fit(tburial1',A_prod2,'gauss1');
    mean_age = f.b1*1e6;
    std_age = f.c1*1e6;

    mc_ages(n) = mean_age;
    mc_errors(n) = std_age;

    if n==(n_runs/100)*floor(n/(n_runs/100))
            progress = 100*n/n_runs;
            clc
            disp(['Progress: ',num2str(progress),'%'])
    end
end % end n_runs

toc

iso_ages = nonzeros(iso_ages)*1e-6;
iso_errors = nonzeros(iso_errors)*1e-6;
mc_ages = nonzeros(mc_ages);
mc_errors = nonzeros(mc_errors);
 
figure(1)
subplot(2,1,1)
histogram(mc_ages*1e-6,edges*1e-6)
hold on
histogram(mc_ages(mc_ages < test_age+dt/2 & mc_ages > test_age-dt/2)*1e-6,edges*1e-6,'FaceColor',[1 0.3 0.3],'FaceAlpha',1)
ax = gca;
% plot([test_age test_age]*1e-6,ax.YLim,'-r','LineWidth',2)
title('Estimated ages (P-PINI)');
xlim([0 3])
xlabel('Age [Myr]')
legend([num2str(n_runs),' runs'],'Correct age')
text(0.05,ax.YLim(2)*0.9,'c)')

subplot(2,1,2)
histogram(iso_ages,edges*1e-6)
hold on
histogram(iso_ages(iso_ages < (test_age+dt/2)*1e-6 & iso_ages > (test_age-dt/2)*1e-6),edges*1e-6,'FaceColor',[1 0.3 0.3],'FaceAlpha',1)
ax = gca;
% plot([test_age test_age]*1e-6,ax.YLim,'-r','LineWidth',2)
title('Estimated ages (Isochron-burial)');
xlim([0 3])
xlabel('Age [Myr]')
legend([num2str(n_runs),' runs'],'Correct age')
text(0.05,ax.YLim(2)*0.9,'d)')

figure(1)
set(gcf,'Units','normalized','position',[0,0.5,0.45,0.4]);

results.mc_ages = mc_ages;
results.iso_ages = iso_ages;
results.mc_errors = mc_errors;
results.iso_errors = iso_errors;
results.iso_MSWD = iso_MSWD;
save('MSWD_test_simple2.mat','results')


figure; plot(iso_ages,iso_MSWD,'.r')
ylabel('MSWD')
xlabel('Age [Myr]')

figure;plot(iso_ages,mc_ages,'*k')
xlabel('Isochron ages')
ylabel('PPINI ages')