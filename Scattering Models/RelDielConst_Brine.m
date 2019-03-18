%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.2: Relative Dielectric Constant of Saline Water
%Description: Code computes the real and imaginary parts of the relative
%dielectric constant of water at any temperature 0<t<30, Salinity 
%0<S<157 0/00, and frequency 0<f<1000GHz
% Using Stogryn (1971) formulations
%Input Variables:
    %t: temperature in C
    %f: frequench in GHz
    %S: salinity in parts per thousand
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%Book Reference: Section 4.2
%MATLAB Code: RelDielConst_SalineWater.m

%Example call: [A B] = RelDielConst_SalineWater(t,f,S)
%Computes the real and imaginary components of the permitivity of saline
    %water based on the temperature value (T) in degrees C, salinity (S) in 0/00
    %and frequency vector (f) and assigns them to vectors A and B respectively

%MATAB CODE

% (C) Jack Landy, 2018 (adapted from code for saline water dielectrics of
% Ulaby et al., 2014)

function [epsr, epsi] = RelDielConst_Brine(T,f,S)
%Note: 'A' matrix used to hold paramaters for each of the equations

%Converts Hz input to GHz as used in this model
f = f*1e9;

% Normality
N = S*(1.707e-2 + 1.205e-5*S + 4.058e-9*S^2);

% Book section 4-1
tau_w = (1.1109e-10 - 3.824e-12*T + 6.938e-14*T^2 - 5.096e-16*T^3)/(2*pi);

eps_w_inf = 4.9;
eps_w_0 = 88.045 - 0.4147*T + 6.295e-4*T^2 + 1.075e-5*T^3;

% Book section 4-5.1
delta = 25-T;

a1 = 1 - 0.255*N + 5.15e-2*N^2 - 6.89e-3*N^3;
b1 = 1 + 0.146e-2*T*N - 4.89e-2*N - 2.97e-2*N^2 + 5.64e-3*N^3;
c1 = 1 - 1.96e-2*delta + 8.08e-5*delta^2 - N*delta*(3.02e-5 + 3.92e-5*delta + N*(1.72e-5 - 6.58e-6*delta));

sigma_b_25 = N*(10.39 - 2.378*N + 0.683*N^2 - 0.135*N^3 + 1.01e-2*N^4);

eps_b0 = eps_w_0*a1;
tau_b = tau_w*b1;
sigma_b = sigma_b_25*c1;

%Complex Permitivity Calculation
eps_0 = 8.854e-12;

epsr = eps_w_inf + (eps_b0 - eps_w_inf)/(1 + (2*pi*f*tau_b)^2);
epsi = (2*pi*f*tau_b)*((eps_b0 - eps_w_inf)/(1 + (2*pi*f*tau_b)^2)) + sigma_b/(2*pi*f*eps_0);

end
