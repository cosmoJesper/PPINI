function [Mean_age,Error,Mean_slope_real,slope_stddev_real,Mean_offset_real,Mean_offset_stddev_real,MSWD] = isochron_mc_LSD_callto(Al,Al_unc,Be,Be_unc,R)
Al = Al*1e-6;
Al(:,2) = Al_unc*1e-6;
Be = Be*1e-6;
Be(:,2) = Be_unc*1e-6;

%ISOCHRON
%Recursive fit to Al-26 and Be-10 data to generate an isochron age.  
%Uses time-varying production rates from Lifton Sato Dunai (2014).  
%Known problems: does not fit line to negative intercepts.  
%Written by Darryl Granger.  Not final--subject to change.  Please ask
%before publishing using this program.


%constants
    %production rate ratio
%     R = 6.97;
%     R = 8.4;
    %Be-10 meanlife
    tau10 = 2.005;
    %Al-26 meanlife
    tau26 = 1.02;
    %effective burial meanlife
    taubur = 1/(1/tau26-1/tau10);
    %SLHL production rate for LSD Be scaling
    P_SLHL = 4.0;
   
    countmax = 20000;
    repmax = 2000;
    
%variables and initial values
    Age = 1;
    hold on;
    rep = 1;
    count = 1;
    Age_diff = 1;
    age_dist = zeros(repmax,1);
    b_dist = zeros(repmax,1);
    m_dist = zeros(repmax,1);
    positives = zeros(repmax,1);
    probability = zeros(repmax,1);
    Be_average = zeros(repmax,1);
    Be_sterr = zeros(repmax,1);
    
%load file AlBedata.dat containing cosmogenic data
%format should be Al-26 Al-26_unc Be-10 Be-10_unc
%units should be in millions of atoms per gram of quartz

% AlBedata = importABfile('Ibe_selected.txt',1,12);
% 
% Al = [AlBedata(:,1), AlBedata(:,2)];
% Be = [AlBedata(:,3), AlBedata(:,4)];

% data = fopen('swiss_foreland2.txt');
% data_all = textscan(data,'%f%f%f%f%s\r\n');
% fclose(data);
% Al = data_all{1}*1e-6;
% Al(:,2) = data_all{2}*1e-6;
% Be = data_all{3}*1e-6;
% Be(:,2) = data_all{4}*1e-6;
% names = data_all{5};
% Al = Al(7:14,:);
% Be = Be(7:14,:);

num_sam = length(Al);

%Create empty arrays for sample-dependent variables
Be_dist = zeros(num_sam,repmax);
Al_dist = zeros(num_sam,repmax);

%Generate local Be-10 production rate and uncertainty using LSD flux
%scaling.  No uncertainty in source elevation--will be added later.

%s_elev = input('what is the uncertainty in your source elevation?  ');

% LSD_DG(lat, lon, elev,0,2000000,1,10);
load LSDout;
P_time = P_SLHL*LSDout.Be;

%Get time series of production rates, smoothed over 100,000 years, and the
%standard deviation of production rates over 101,000 years.  
%Offset time by 100,000 years for exposure on the landscape prior to burial.

P_smooth = smooth(P_time,101);
P_std = movingstd(P_time,50,'c');
time = LSDout.tv - 50000;

%Make an initial production rate estimate, using the long-term mean.
P = mean(P_smooth);

%Generate normalized Al and Be data
Al_star = [Al(:,1)./(P*tau26*R), Al(:,2)./(P*tau26*R)];
Be_star = [Be(:,1)./(P*tau10), Be(:,2)./(P*tau10)];

%Regress a slope using York's method
    X = Be_star(:,1);
    Y = Al_star(:,1);
    %X = AlBedata(:,3);
    %Y = AlBedata(:,1);

    sX = Be_star(:,2);
    sY = Al_star(:,2);
    
    %sX = AlBedata(:,4);
    %sY = AlBedata(:,2);

    [b,m,sb,sm,S]  =york_fit(X',Y',sX',sY');
    N=numel(X);

 %estimate initial age from slope
    Age = -taubur*log(m/2);

    for rep = 1:repmax
        
        %Randomize data by [production rate (boxcar)--later]  measurement errors
        %(Gaussian).
        
        Al_mc = Al(:,1) + randn*Al(:,2);
        Be_mc = Be(:,1) + randn*Be(:,2);
        P_mc = randn;
        
        
       for count=1:countmax

           %POSTBURIAL 
           %Calculate postburial production rate ratio assuming constant
            %production for entire burial time (choose this or next line)
%             R_post = (1-exp(-Age/tau26))/(1-exp(-Age/tau10));
            R_post = R;
            
            %Calculate postburial production rate ratio for recent exposure
            %(choose this or previous line)
%             R_post = 2;
%             R_post = 6.97;
            
            %Find intercept of isochron with postburial ratio line
            Be_post = b/(R_post - m);
            Al_post = Be_post*R_post;

           if or(Be_post<0,Al_post<0)
               Be_post = 0;
               Al_post = 0;
           end;

            %Subtract postburial production from Al and Be data to get inherited
            %component.  If calculated postburial is greater than measured values, 
            %then the postburial component is limited to the lowest
            %measurement.
            Al_inh =  Al_star(:,1) - Al_post;
            Be_inh =  Be_star(:,1) - Be_post;

            if (min(Be_inh)<0)
            Be_post = min(Be_star(:,1));
            Al_post = Be_post*R_post;
            Al_inh =   Al_star(:,1) - Al_post;
            Be_inh =  Be_star(:,1) - Be_post;
            end

            if (min(Al_inh)<0)
            Al_post = min(Al_star(:,1));
            Be_post = Al_post/R_post;
            Al_inh = Al_star(:,1) - Al_post;
            Be_inh = Be_star(:,1) - Be_post;
            end
        
        %PREBURIAL
        %Generate Be-10 linearization factor based on erosion rates at time zero
        %Use first line for constant exposure OR second line for steady
        %erosion.
%         Be_corr = 1./(1+Be_inh.*exp(Age/tau10));
%         Be_corr = 1;
        Be_corr = 1 - 0.5.*Be_inh*exp(Age/tau10);
        Be_lin = Be_inh.*Be_corr;

        %Add back in the postburial production for fitting
        Be_fit = Be_lin + Be_post;
    
        %Make placeholder for earlier age
        Old_age = Age;

        %Regress a slope using York's method.  
         X = Be_fit;
         Y = Al_star(:,1);
        sX = Be_star(:,2);
        sY = Al_star(:,2);
       [b,m,sb,sm,S]  = york_fit(X',Y',sX',sY');

        %estimate age from slope
        Age = -taubur*log(abs(m)*tau26/tau10);
        
        %Get production rate at time zero using new age
        P_age = interp1(time,P_time,real(Age));
        sP = interp1(time, P_std, Age);
%         P_age = P(1);
        
        %randomize the new production rate, using the Gaussian random
        %factor generated at the beginning of the outer loop.
        P = P_age + P_mc*P_std;
       
        %Generate new normalized Al and Be data
         Al_star = [Al_mc(:,1)./(P*tau26*R), Al(:,2)./(P*tau26*R)];
         Be_star = [Be_mc(:,1)./(P*tau10), Be(:,2)./(P*tau10)];
         
         Age_diff = abs(Old_age - Age);
         
        if Age_diff<0.005
            break
        end;
    end;
    
   if count>countmax-1
       'AGE DID NOT CONVERGE'
   end
    
    if b>0
        positives(rep) = 1;
    end;
    
   Be_dist(:,rep) = Be_fit;
   Al_dist(:,rep) = Al_star(:,1);
   Be_post_dist(:,rep) = Be_post;
   Al_post_dist(:,rep) = Al_post;
   P_dist(:,rep) = P;
    
   m_dist(rep) = m;
   b_dist(rep) = b;
   age_dist(rep) = Age;
   S_dist(rep) = S;
   probability(rep) = chi2pdf(S,(N-2));
   
   %Get average 10Be concentration at time zero
   Be_original = Be_inh*exp(Age/tau10);
   Be_average(rep) = mean(Be_original);
   Be_sterr(rep) = std(Be_original)/sqrt(num_sam);
   
    end;


%plot error ellipse for each point

%for i=1:N
 %   [Xe Ye] = ellipse(mean(Be_dist(i,:)),mean(Al_dist(i,:)),std(Be_dist(i,:)),std(Al_dist(i,:)),32);
  %  plot(Xe,Ye,'k')
%end

%Generate curves from fit
%normalized values reported
Mean_Be_post = mean(Be_post_dist);
Mean_Al_post = mean(Al_post_dist);
s_Be_post = std(Be_post_dist);
s_Al_post = std(Al_post_dist);
Mean_slope = mean(m_dist);
slope_error = std(m_dist);
Mean_offset = mean(b_dist);
offset_error = std(b_dist);
Mean_age = mean(age_dist);
Positives_age = sum(age_dist.*positives)/sum(positives);
Error = std(age_dist);
MSWD = mean(S_dist)/(num_sam-2);
Probability = mean(probability);


%not-normalized values reported
Mean_Be_post_real = Mean_Be_post*P*tau10*10^6;
Mean_Al_post_real = Mean_Al_post*P*tau26*R*10^6;
Mean_offset_real = mean(Mean_offset*P*tau26*R*10^6);
Mean_offset_stddev_real = mean(offset_error*P*tau26*R*10^6);
Mean_slope_real = mean(m_dist)*P*tau26*R/(P*tau10);
slope_stddev_real = std(m_dist)*P*tau26*R/(P*tau10);

Mean_Be_post_real = Mean_Be_post*P*tau10*10^6;
Mean_Be_post_stddev_real = s_Be_post*P*tau10*10^6;

Mean_age = mean(age_dist);
Error = std(age_dist);
MSWD = mean(S_dist)/(num_sam-2);
Probability = mean(probability);

% 
% %x1 = Be_star(:,1)*(P*tau10);
% %y1 = Al_star(:,1)*(P*tau26*R);
% %sx1 = Be_star(:,2)*(P*tau10);
% %sy1 = Al_star(:,2)*(P*tau26*R);
% 
% P*tau26*R/(P*tau10)
% 
% %write file
% fileID = fopen('results.txt', 'w');
% fprintf(fileID, 'Mean_Be_post\t %e\n', Mean_Be_post);
% fprintf(fileID, 's_Be_post\t %e\n', s_Be_post);
% fprintf(fileID, 'Mean_Al_post\t %e\n', Mean_Al_post);
% fprintf(fileID, 's_Al_post\t %e\n', s_Al_post);
% fprintf(fileID, 'Mean_slope\t %e\n', Mean_slope);
% fprintf(fileID, 'slope_error\t %e\n', slope_error);
% fprintf(fileID, 'Mean_offset\t %e\n', Mean_offset);
% fprintf(fileID, 'offset_error\t %e\n', offset_error);
% fprintf(fileID, 'Mean_age\t %e\n', Mean_age);
% fprintf(fileID, 'Positives_age\t %e\n', Positives_age);
% fprintf(fileID, 'Error\t %e\n', Error);
% fprintf(fileID, 'MSWD\t %e\n', MSWD);
% fprintf(fileID, 'Probability\t %e\n\n\n', Probability);
% fprintf(fileID, 'Values NOT normalized\t \n');
% fprintf(fileID, 'linearization factor\t %e\n', Be_corr);
% 
% 
% fprintf(fileID, 'Mean_Be_post_real\t %e\n', Mean_Be_post_real);
% fprintf(fileID, 'Mean_Al_post_real\t %e\n', Mean_Al_post_real);
% fprintf(fileID, 'Mean_offset_real\t %e\n', Mean_offset_real);
% fprintf(fileID, 'Mean_offset_stddev_real\t %e\n', Mean_offset_stddev_real);
% fprintf(fileID, 'Mean_slope_real\t %e\n', Mean_slope_real);
% fprintf(fileID, 'slope_stddev_real\t %e\n', slope_stddev_real);
% fprintf(fileID, 'Mean_Be_post_real\t %e\n', Mean_Be_post_real);
% fprintf(fileID, 'Mean_Be_post_stddev_real\t %e\n', Mean_Be_post_stddev_real);
% fprintf(fileID, 'Mean_age\t %e\n', Mean_age);
% fprintf(fileID, 'Error\t %e\n', Error);
% fclose(fileID);
% 
% 
% Al_star = [Al(:,1)./(P*tau26*R), Al(:,2)./(P*tau26*R)];
% Be_star = [Be(:,1)./(P*tau10), Be(:,2)./(P*tau10)];
% 
% xx1 = Be_star(:,1);
% yy1 = Al_star(:,1);
% sxx1 = Be_star(:,2);
% syy1 = Al_star(:,2);
% 
% 
% for i=1:N
%     [Xe Ye] = ellipse(xx1(i),yy1(i),sxx1(i),syy1(i),32);
%     plot(Xe,Ye,'r')
% 
%     [Xe Ye] = ellipse(xx1(i)*Be_corr(i),yy1(i),sxx1(i),syy1(i),32);
%     plot(Xe,Ye,'k')
% end
% xlabel('Be-10*');
% ylabel('Al-26*');
% title('Isochron diagram');
% 
% 
% x2 = 0:0.0001:(max(xx1)+0.1*max(xx1));
% y_best= Mean_slope*x2+Mean_offset;
% y_high= (Mean_slope+slope_error)*x2+Mean_offset+offset_error;
% y_low= (Mean_slope-slope_error)*x2+Mean_offset-offset_error;
% 
% %y_best =2*exp(-Mean_age/taubur)*(x2-Mean_Be_post)./(1 + x2*exp(Mean_age/taubur))+Mean_Al_post;
% %y_high =2*exp(-(Mean_age-Error)/taubur)*(x2-(Mean_Be_post + offset_error/R_post))./(1 + x2*exp((Mean_age-Error)/taubur))+(Mean_Al_post + offset_error);
% %y_low =2*exp(-(Mean_age+Error)/taubur)*(x2-(Mean_Be_post - offset_error/R_post))./(1 + x2*exp((Mean_age+Error)/taubur))+(Mean_Al_post - offset_error);
% 
% %plot(mean(Be_dist'),mean(Al_dist'),'k+')
% 
% plot(x2,y_best)
% plot(x2,y_high)
% plot(x2,y_low)
% 
% hold off;



