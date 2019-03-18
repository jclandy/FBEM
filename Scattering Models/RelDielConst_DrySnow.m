%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.5: Relative Dielectric Constant of Dry Snow
%Description: Code computes the real and imaginary parts of the relative
%dielectric constant of Dry Snow
%Input Variables:
    %T: Temperature in C
    %Ps: Dry Snow Density in g/cm^3
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%Book Reference: Section 4-6
%MATLAB Code: RelDielConst_DrySnow.m

%Example call: [A B] = RelDielConst_DrySnow(T,Ps,f)
%Computes the real and imaginary components of the permitivity of Dry Snow
    %based on the temperature value (T) in degrees C, dry snow density (Ps), and frequency
    %vector (f) and assigns them to vectors A and B respectively

%MATAB CODE

function [epsr_ds,epsi_ds] = RelDielConst_DrySnow(T, Ps, f )

vi=  Ps /0.9167 ;


[epsr_ice, epsi_ice] = RelDielConst_PureIce(T,f);

if vi <= 0.45
    epsr_ds = 1 + 1.4667 .* vi + 1.435 .* vi.^3; % 0<vi<0.45
end
if vi > 0.45
    epsr_ds = (1+ 0.4759 .*vi).^3; % vi>0.45
end

epsi_ds = 0.34 * vi * epsi_ice ./(1- 0.42 * vi).^2;


end


%Code 4.3: Relative Dielectric Constant of Pure Ice
%Description: Code computes the real and imaginary parts of the relative
%dielectric constant of pure ice at any temperature -40C<T<0C and 
%frequency 1<f<1000GHz
%Input Variables:
    %T: Temperature in degree C
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant
    %epsi: imaginary part of relative dielectric constant
%Book Reference: Section 4-3
%MATLAB Code: RelDielConst_PureIce.m

%Example call: [A B] = RelDielConst_PureIce(T,f)
%Computes the real and imaginary components of the permitivity of pure ice
    %based on the temperature value (T) in degrees C and frequency
    %vector (f) and assigns them to vectors A and B respectively

%MATLAB CODE

function [epsr epsi] = RelDielConst_PureIce(T,f)

T = T + 273; % represent temperature in Kelvin

theta = (300/T) - 1;

B1 = 0.0207;
B2 = 1.16e-11;
b = 335;

alpha = (0.00504 + 0.0062*theta)*exp(-22.1*theta);

betaM = (B1/T) * exp(b/T) / (exp(b/T)-1)^2  + B2.*f.^2;
delBeta = exp(-9.963 + 0.0372*(T-273.16));

beta = betaM + delBeta;

% epsr = 3.1884 + 9.1e-4*(T-273) .*f./f;

epsr = 3.1884 + 9.1e-4 *(T-273);

epsi = alpha./f + beta.*f;
end
