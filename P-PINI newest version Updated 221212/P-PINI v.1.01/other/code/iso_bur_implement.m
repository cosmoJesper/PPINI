function [age_final,error,mean_slope_real,mean_offset_real] = iso_bur_implement(Al,Al_unc,Be,Be_unc,PAl,PBe)

NAl = Al;
NAl_unc = Al_unc;
NBe = Be;
NBe_unc = Be_unc;
P_Al = PAl;
P_Be = PBe;

if numel(NAl_unc) == 1
    NAl_unc = NAl_unc*ones(size(NAl));
end

if numel(NBe_unc) == 1
    NBe_unc = NBe_unc*ones(size(NBe));
end

%%%%%%%%%%%%%%%%%%%%
n = numel(NAl);
m_dist = zeros(1,n);
b_dist = zeros(1,n);
sigma_age_dist = zeros(1,n);
age_dist = zeros(1,n);
%%%%%%%%%%%%%%%%%%%%

R_post = 6.97;
% R_post = 2;

TBe = 1.387e6;
TAl = 0.705e6;
lambda_Be = log(2)/TBe;
lambda_Al = log(2)/TAl;
tauBe = 1/lambda_Be;
tauAl = 1/lambda_Al;

Al_norm = NAl/(P_Al*tauAl);
Al_unc_norm = NAl_unc/(P_Al*tauAl);
Be_norm = NBe/(P_Be*tauBe);
Be_unc_norm = NBe_unc/(P_Be*tauBe);

Nmax = 200;
tol = 10e-18;
m = 0.005;

w_Al = 1./Al_unc_norm.^2;
w_Be = 1./Be_unc_norm.^2;

ms = zeros(1,Nmax+1);
ms(1) = m;

for i = 1:Nmax
    W = w_Al.*w_Be./(m^2*w_Al+w_Be);
    
    Be_bar = sum(W.*Be_norm)/sum(W);
    Al_bar = sum(W.*Al_norm)/sum(W);
    
    u = Be_norm - Be_bar;
    v = Al_norm - Al_bar;
    
    beta = W.*(u./w_Al + m*v./w_Be);
    m = sum(W.*beta.*v)/sum(W.*beta.*u);
    
    ms(i+1) = m;
    if abs((ms(i+1)-ms(i))/ms(i+1)) < tol
        break
    end
end

b = Al_bar - m*Be_bar;


x = Be_bar+beta;
x_bar = sum(W.*x)/sum(W);
u = x - x_bar;
% dev = Al_norm - m*Be_norm - b;


sigma_m = sqrt(1/sum(W.*u.^2));
% sigma_b = sqrt(1./sum(W)+x_bar^2*sigma_m^2);


age = -log(m/2)/(lambda_Al - lambda_Be);

for rep = 1:1000
    Al_mc = NAl + randn*NAl_unc;
    Be_mc = NBe + randn*NBe_unc;
    P_Be_mc = P_Be + randn*0.05;
    P_Al_mc = P_Al + randn*0.05;

    for j = 1:50000

        Be_post = b/(R_post - m);
        Al_post = Be_post*R_post;

        if Be_post < 0
            Be_post = 0;
            Al_post = 0;
        end

        Al_inh = Al_norm - Al_post;
        Be_inh = Be_norm - Be_post;

        if min(Be_inh < 0)
            Be_post = min(Be_norm);
            Al_post = Be_post*R_post;
            Al_inh = Al_inh - Al_post;
            Be_inh = Be_inh - Be_post;
        end

        if min(Al_inh < 0)
            Al_post = min(Al_norm);
            Be_post = Al_post/R_post;
            Al_inh = Al_inh - Al_post;
            Be_inh = Be_inh - Be_post;
        end    
        
        Be_corr = 1;
%         Be_corr = 1- 0.5*Be_inh*exp(age/tauBe);
        Be_lin = Be_inh.*Be_corr;

        Be_fit = Be_lin + Be_post;

        age_old = age;


    %%%%%% Doing york %%%%%%%%%
        w_Al = 1./Al_unc_norm.^2;
        w_Be = 1./Be_unc_norm.^2;

        ms = zeros(1,Nmax+1);
        ms(1) = m;

        for i = 1:Nmax
            W = w_Al.*w_Be./(m^2*w_Al+w_Be);

            Be_bar = sum(W.*Be_fit)/sum(W);
            Al_bar = sum(W.*Al_norm)/sum(W);

            u = Be_fit - Be_bar;
            v = Al_norm - Al_bar;

            beta = W.*(u./w_Al + m*v./w_Be);
            m = sum(W.*beta.*v)/sum(W.*beta.*u);

            ms(i+1) = m;
            if abs((ms(i+1)-ms(i))/ms(i+1)) < tol
                break
            end
        end

        b = Al_bar - m*Be_bar;


        x = Be_bar+beta;
        x_bar = sum(W.*x)/sum(W);
        u = x - x_bar;
        dev = Al_norm - m*Be_fit - b;
        S = sum(W.*dev.^2);


        sigma_m = sqrt(1/sum(W.*u.^2));
%         sigma_b = sqrt(1./sum(W)+x_bar^2*sigma_m^2);
    %%% not doing york anymore %%%%%

        age = -log(m*(tauAl/tauBe))/(lambda_Al-lambda_Be);
        age_max = -log((m-sigma_m)*(tauAl/tauBe))/(lambda_Al-lambda_Be);
        age_min = -log((m+sigma_m)*(tauAl/tauBe))/(lambda_Al-lambda_Be);

        Al_norm = Al_mc/(P_Al_mc*tauAl);
        Al_unc_norm = NAl_unc/(P_Al_mc*tauAl);
        Be_norm = Be_mc/(P_Be_mc*tauBe);
        Be_unc_norm = NBe_unc/(P_Be_mc*tauBe);

        if abs(age_old - age) < 0.005
            break
        end
    end
    
    if j>49999
        'Not converging'
    end
        
   m_dist(rep) = m;
   b_dist(rep) = b;
   sigma_age_dist(rep) = (age_max - age_min)/2;
   age_dist(rep) = age;
end

mean_offset_real = mean(b_dist)*P_Al*tauAl;
mean_slope_real = mean(m_dist)*P_Al*tauAl/(P_Be*tauBe);

age_final = mean(age_dist);
error = std(age_dist);
% error = sigma_age_dist(rep);


