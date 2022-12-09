function out = LSDscaling(h,Rc,SPhi,w,consts,nuclide)

% Implements the Lifton Sato Dunai scaling scheme for spallation.
%
% Syntax: scalingfactor = LiftonSatoSX(h,Rc,SPhi,w,consts);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%   
%
% Vectorized. Send in scalars or vectors of common length. 
%

% Written by Nat Lifton 2013, Purdue University
% nlifton@purdue.edu
% Based on code by Greg Balco -- Berkeley Geochronology Lab
% balcs@bgc.org
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2001-2013, University of Washington, Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).


% what version is this?

ver = '1.0';

mfluxRef = consts.mfluxRef;
muRef = (mfluxRef.neg + mfluxRef.pos);

% Select reference values for nuclide of interest or flux

if nuclide == 3
    HeRef = consts.P3nRef + consts.P3pRef;
elseif nuclide == 10
    BeRef = consts.P10nRef + consts.P10pRef;
elseif nuclide == 14
    CRef = consts.P14nRef + consts.P14pRef;
elseif nuclide == 26
    AlRef = consts.P26nRef + consts.P26pRef;
else    
    SpRef = consts.nfluxRef + consts.pfluxRef;
    % Sato et al. (2008) Reference hadron flux integral >1 MeV
end

EthRef = consts.ethfluxRef;
ThRef = consts.thfluxRef;

% Site nucleon fluxes

NSite = Neutrons(h,Rc,SPhi,w,consts,nuclide);
[ethflux,thflux] = NeutronsLowE(h,Rc,SPhi,w);
PSite = Protons(h,Rc,SPhi,consts,nuclide);

% Site omnidirectional muon flux
mflux = Muons(h,Rc,SPhi);%Generates muon flux at site from Sato et al. (2008) model
muSite = (mflux.neg + mflux.pos);

%Nuclide-specific scaling factors as f(Rc)
if nuclide == 3
    Site.He = (NSite.P3n + PSite.P3p)./HeRef;
elseif nuclide == 10
    Site.Be = (NSite.P10n + PSite.P10p)./BeRef;
elseif nuclide == 14
    Site.C = (NSite.P14n + PSite.P14p)./CRef;
elseif nuclide == 26
    Site.Al = (NSite.P26n + PSite.P26p)./AlRef;
else    %Total nucleon flux scaling factors as f(Rc)
    Site.sp = ((NSite.nflux + PSite.pflux))./SpRef; % Sato et al. (2008) Reference hadron flux integral >1 MeV
end

Site.E = NSite.E;%Nucleon flux energy bins
Site.eth = ethflux./EthRef; %Epithermal neutron flux scaling factor as f(Rc)
Site.th = thflux./ThRef;%Thermal neutron flux scaling factor as f(Rc)

%Differential muon flux scaling factors as f(Energy, Rc)
Site.muE = mflux.E;%Muon flux energy bins (in MeV)
Site.mup = mflux.p;%Muon flux momentum bins (in MeV/c)

for i = 1:length(Rc)
    Site.muSF(i,:) = muSite(i,:)./muRef;
end

%Integral muon flux scaling factors as f(Rc)
Site.muTotal = mflux.total./mfluxRef.total;%Integral total muon flux scaling factor
Site.mn = mflux.nint./mfluxRef.nint;%Integral neg muon flux scaling factor
Site.mp = mflux.pint./mfluxRef.pint;%Integral pos muon flux scaling factor
Site.mnabs = mflux.nint;%Integral neg muon flux
Site.mpabs = mflux.pint;%Integral pos muon flux 
     
out = Site;

