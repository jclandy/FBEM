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
%0<S<40 0/00, and frequency 0<f<1000GHz
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

function [epsr epsi] = RelDielConst_SalineWater(t,f,S)
%Note: 'A' matrix used to hold paramaters for each of the equations

%Converts Hz input to GHz as used in this model

%Conductvity
A = [2.903602 8.607e-2 4.738817e-4 -2.991e-6 4.3041e-9];
sig35 = A(1) + A(2)*t + A(3)*t^2 + A(4)*t^3 + A(5)*t^4;

A = [37.5109 5.45216 0.014409 1004.75 182.283];
P = S * ((A(1) + A(2)*S + A(3)*S^2) / (A(4) + A(5)*S + S^2));

A = [6.9431 3.2841 -0.099486 84.85 69.024];
alpha0 = (A(1) + A(2)*S + A(3)*S^2) / (A(4) + A(5)*S + S^2);

A = [49.843 -0.2276 0.00198];
alpha1 = A(1) + A(2)*S + A(3)*S^2;

Q = 1 + ((alpha0*(t-15))/(t+alpha1));

sigma = sig35*P*Q;


%Other Model Paramaters
a=[0.46606917e-2 -0.26087876e-4 -0.63926782e-5 0.63000075e1 0.26242021e-2 -0.42984155e-2 ...
   0.34414691e-4 0.17667420e-3 -0.20491560e-6 0.58366888e3 0.12634992e3 0.69227972e-4 ...
   0.38957681e-6 0.30742330e3 0.12634992e3 0.37245044e1 0.92609781e-2 -0.26093754e-1];

epsS = 87.85306*exp(-0.00456992*t - a(1)*S - a(2)*S^2 - a(3)*S*t);
epsOne = a(4)*exp(-a(5)*t-a(6)*S-a(7)*S*t);
tau1 = (a(8)+a(9)*S)*exp(a(10)/(t+a(11)));
tau2 = (a(12)+a(13)*S)*exp(a(14)/(t+a(15)));
epsInf = a(16) + a(17)*t + a(18)*S;

%Complex Permitivity Calculation
eps = ((epsS-epsOne)./(1-1i*2*pi.*f.*tau1)) + ((epsOne-epsInf)./(1-1i*2*pi.*f.*tau2)) + (epsInf) + 1i*((17.9751*sigma)./f);

%Seperates real and imaginary components
epsr = real(eps);
epsi = imag(eps);

end
