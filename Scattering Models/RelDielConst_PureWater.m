%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.1: Relative Dielectric Constant of Pure Water
%Description: Code computes the real and imaginary parts of the relative
    %dielectric constant of water at any temperature 0<t<30 and frequency
    %0<f<1000 GHz. Uses the double-Debye model.
%Input Variables:
    %t: temperature in C
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%Book Reference: Section 4.1
%MATLAB Code:  RelDielConst_PureWater.m

%Example call: [A B] = RelDielConst_PureWater(t,f)
%Computes the real and imaginary components of the permitivity of pure
    %water based on the temperature value (t) in degrees C and frequency
    %vector (f) and assigns them to vectors A and B respectively

%MATLAB Code

function [epsr,epsi] = RelDielConst_PureWater(t,f)
   
    %Model Paramaters
        a=[0.63000075e1 0.26242021e-2 0.17667420e-3 0.58366888e3 0.12634992e3 ...
            0.69227972e-4 0.30742330e3 0.12634992e3 0.37245044e1 0.92609781e-2]; %Eqn parameters

        epsS = 87.85306*exp(-0.00456992*t);
        epsOne = a(1)*exp(-a(2)*t);
        tau1 = a(3)*exp(a(4)/(t+a(5)));
        tau2 = a(6)*exp(a(7)/(t+a(8)));
        epsInf = a(9) + a(10)*t;

        %Permitivity Calculation
        eps = ((epsS-epsOne)./(1-1j*2*pi.*f.*tau1)) + ((epsOne-epsInf)./(1-1j*2*pi.*f.*tau2)) + (epsInf);

        %Extracts the real and imaginary components of relative permitivity
        epsr = real(eps);
        epsi = imag(eps);
  end