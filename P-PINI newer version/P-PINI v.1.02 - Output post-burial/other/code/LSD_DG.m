function LSD_DG(lat,lon,alt,atm,age,w,nuclide)

% This function calculates Lifton, Sato, and Dunai time-dependent scaling factors 
% for a given set of inputs
% syntax : LSD(lat,lon,alt,atm,age,nuclide);

% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes, 
%     or 0-360 degrees E) 
% alt = sample altitude in m above sea level
% atm = atmospheric model to use: 1 for U.S. Standard Atmosphere, 
%     0 for ERA-40 Reanalysis
% age = age of sample
% w = gravimetric fractional water content - 0.066 is default 
%     typically about 14% volumetric per Fred Phillips. -1 gives default
%     value
% nuclide = nuclide of interest: 26 for 26Al, 10 for 10Be, 14 for 14C,
%     3 for 3He, 0 for nucleon flux
% 
% Input values as scalars
%
% Based on code written by Greg Balco -- Berkeley
% Geochronology Center
% balcs@bgc.org
% 
% Modified by Brent Goehring and 
% Nat Lifton -- Purdue University
% nlifton@purdue.edu, bgoehrin@purdue.edu
%
% Copyright 2013, Berkeley Geochronology Center and
% Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).

% Modified by Darryl Granger to simplify for burial and paleo-erosion rate
% calculations.  Changed to uniform 1000-year spacing.

% what version is this?
ver = '1.0_DG';

load consts_LSD;

is14 = 0;
is10 = 0;
is26 = 0;
is3 = 0;
isflux = 0;

% Load the input data structure

sample.lat = lat;
sample.lon = lon;
sample.alt = alt;
sample.atm = atm;
sample.age = age;
sample.nuclide = nuclide;

if nuclide == 14
    is14 = 1;
elseif nuclide == 10
    is10 = 1;
elseif nuclide == 26
    is26 = 1;    
elseif nuclide == 3
    is3 = 1;  
else
    isflux = 1;
end    

if sample.atm == 1
    stdatm = 1;
    gmr = -0.03417; % Assorted constants
    dtdz = 0.0065; % Lapse rate from standard atmosphere
else
    stdatm = 0;
end
    
% Make the time vector
calFlag = 0;

% Age Relative to t0=2010
%DG: changed to uniform 1000 year spacing (adding 60 to start in 1950)
tv = [60:1000:2000060 logspace(log10(2001060),7,200)];

LSDRc = zeros(1,length(tv));

% Need solar modulation parameter
%DG: changed to make SPhi constant.
this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
%this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

if w < 0
    w = 0.066; % default gravimetric water content for Sato et al. (2008)
end

% interpolate an M for tv > 7000...
%DG: changed to interpolate over whole range of tv
temp_M = interp1(consts.t_M,consts.M,tv);

% Pressure correction

if stdatm == 1
    % Calculate site pressure using the Standard Atmosphere parameters with the
    % standard atmosphere equation.
    sample.pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (alt.*dtdz)) ) );
else    
    sample.pressure = ERA40atm(sample.lat,sample.lon,sample.alt);
end

% catch for negative longitudes before Rc interpolation
if sample.lon < 0; sample.lon = sample.lon + 360;end;

% Make up the Rc vectors.

% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average. 

dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

LSDRc = temp_M.*(dd(1)*cosd(sample.lat) + ...
   dd(2)*(cosd(sample.lat)).^2 + ...
   dd(3)*(cosd(sample.lat)).^3 + ...
   dd(4)*(cosd(sample.lat)).^4 + ...
   dd(5)*(cosd(sample.lat)).^5 + ...
   dd(6)*(cosd(sample.lat)).^6); 
% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
%DG: change interpolation range and reverse order of this and previous section

[loni,lati,tvi] = meshgrid(sample.lon,sample.lat,tv(1:7));
LSDRc(1:7) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,loni,lati,tvi);

% Next, chop off tv
%DG: calculate whole range.  comment out this section
%clipindex = find(tv <= sample.age, 1, 'last' );
%tv2 = tv(1:clipindex);
%if tv2(end) < sample.age;
%    tv2 = [tv2 sample.age];
%end;
% Now shorten the Rc's commensurately 
%LSDRc = interp1(tv,LSDRc,tv2);
%DG: change to match tv
LSDSPhi = this_SPhi;

LSDout = LSDscaling(sample.pressure,LSDRc(:),LSDSPhi,w,consts,nuclide);

%DG: change to tv
LSDout.tv = tv;
LSDout.Rc = LSDRc;
LSDout.pressure = sample.pressure;
LSDout.alt = sample.alt;

% Write results to file

save ('LSDout', 'LSDout');

%disp('Output will be saved in native binary Matlab format as a .mat file');
%name = input('Give the output a unique name to identify it: ','s');
%save(name,'LSDout');

%CD = pwd;
%disp(['Output saved to ' CD]) 

%clear all; close all;

