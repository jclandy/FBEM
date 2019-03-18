%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.4: TVB Dielectric Model for Heterogeneous Mixtures

%Description: Code computes the Tinga-Voss-Blossey (TVB) model for a 
% heterogeneous mixture composed of inclusions with eps_i in a host material 
% with eps_h. The inclusions are randomly oriented, but their shapes can be 
% specified. 
    
%Input Variables:
    %eps_i: complex dielectric constant of inclusion material
    %eps_h: complex dielectric constant of host material
    %shape: shape of the inclusion
        % 1: circular disc
        % 2: spherical 
        % 3: needle
    % vi: inclusion volume fraction
    
%Output Products:
    %eps_m: complex dielectric constant of mixture

%Book Reference: Section 4-4.3

%Example call: [A ] = TVBmodel_HeterogeneousMix(eps_i, eps_h, shape, vi)
    
%MATLAB Code

function [eps_m] = TVBmodel_HeterogeneousMix(eps_i, eps_h, shape, vi)

if shape == 1 %case of thin circular disc inclusions
    eps_m = eps_h + vi/3*(eps_i - eps_h)*(2*eps_i*(1-vi)+eps_h*(1+2*vi))...
        ./(vi*eps_h +(1-vi)*eps_i);
end

if shape ==2 % case of spherical inclusions
    eps_m = eps_h + 3*vi*eps_h*(eps_i -eps_h)./((2*eps_h+eps_i)-vi*(eps_i-eps_h));
end

if shape == 3 % case of needle inclusions
    eps_m = eps_h + vi/3*(eps_i-eps_h)*(eps_h*(5+vi)+(1-vi)*eps_i)./ ...
        (eps_h*(1+ vi)+eps_i*(1-vi));
end

  end