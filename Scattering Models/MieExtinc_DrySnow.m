%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 11.3: Mie Extinction of Dry Snow

%Description: Code computes the Mie absorption, scattering, and extinction
%coefficients of dry snow for specified snow density and ice particle
%radius.

%Input Variables:
    % rho_s: Snow density (g/cm3)
    % ri: Ice-particle radius (m)
    % f: frequency (GHz)
    % t: temperature (degree C)
    
%Output Products:
    % kappa_a: absorption coefficient (Np/m)
    % kappa_s: scattering coefficient (Np/m)
    % kappa_e: extinction coefficient (Np/m)
    % a: single-scattering albedo
    
%Book Reference: Section 11-15.1 and eq 11.112 and 11.113 with Q computed
%according to the Mie model of section 8-5.

%Matlab code

function [kappa_a,kappa_s,kappa_e,a] = MieExtinc_DrySnow(rho_s,ri,f,t)

eps_b = 1; % dielectric constant of background
rho_i = 0.9167; %density of ice (g /cm3)

%- calculate relative dielectric constant of pure ice
[epsr epsi] = RelDielConst_PureIce(t,f);
eps_sp = epsr -1i * epsi;

%-- calculate the Mie efficiencies using Code 8.12
[Es Ea Ee Eb t1 t2 t3 t4] = ...
    Mie_Rayleigh_ScatteringOfSpheres(ri, f, eps_sp, eps_b); 

area = pi * ri^2;

Qs = Es * area;
Qa = Ea * area;
Qe = Ee * area;

Nv = rho_s ./rho_i ./(4/3 * pi .*ri^3); 

kappa_a = Nv .* Qa;
kappa_s = Nv .* Qs; 
kappa_e = Nv .* Qe;

a = kappa_s /kappa_e;

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
    

%Code 8.12: Mie and Rayleigh Scattering By Spherical Particle

%Description: Code computes the absorption, scattering, extinction, and
%backscattering efficiencies of a dielectric sphere according to bot the
%Mie solution and the Rayleigh approximation.

%Input Variables
    %r: radius of particle (meters)
    %f: frequency (GHz)
    %epsp: epsp1-j*epsp2 for particles
    %epsb = epsb1 -j epsb2 for background medium 
    
%Output Variables
    %Es: Mie Scattering Efficiency 
    %Ea: Mie Absorption Efficiency
    %Ee: Mie Extinction Efficiency
    %Eb: Mie Backscattering Efficiency
    
    %eta_s_r: Rayleigh Scattering Efficiency 
    %eta_a_r: Rayleigh Absorption Efficiency
    %eta_e_r: Rayleigh Extinction Efficiency
    %eta_b_r: Rayleigh Backscattering Efficiency

%Matlab Code

function [Es Ea Ee Eb eta_s_r eta_a_r eta_e_r eta_b_r] = ... 
    Mie_Rayleigh_ScatteringOfSpheres(r, f, epsp, epsb)

epsb1 = real(epsb);

np = sqrt(epsp); % index of refraction of spherical particle
nb = sqrt(epsb); % index of refraction of background medium

n = np ./nb; % relative index of refraction
%  n = np ./sqrt(epsb1);


chi = 20./3 *pi*r*f*sqrt(epsb1); % normalized circumference in background

%--- Calculate Rayleigh Approximation solution

BigK = (n.^2 - 1) ./(n.^2 + 2); 

eta_s_r = 8/3 .* chi.^4 .* abs(BigK).^2;

eta_a_r = 4 .* chi .* imag(-BigK); 

eta_e_r = eta_a_r + eta_s_r; % Extinction efficiency

eta_b_r = 4 * chi.^4 .* abs(BigK).^2;  % backscattering efficiency



%--- Calculate Mie Scattering solution

%Calculation of Es
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            
    %A_1=    (sin(real(n).*chi).*cos(real(n).*chi)+j*sinh(imag(n).*chi).*cosh(imag(n).*chi))...
    %        ./(sin(real(n).*chi).^2 + sinh(imag(n).*chi).^2);

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (2*l + 1).*(abs(a).^2+abs(b).^2);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end
Es = 2./(chi).^2 .* runSum;

%Calculation of Ee
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            
    %A_1 =   (sin(real(n).*chi).*cos(real(n).*chi)+j*sinh(imag(n).*chi).*cosh(imag(n).*chi))...
    %       ./(sin(real(n).*chi).^2 + sinh(imag(n).*chi).^2);

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (2*l + 1).*real(a+b);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end
Ee = 2./(chi).^2 .* runSum;

%Calculation of Eb
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (-1)^l .* (2*l + 1) .*(a-b);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end

Eb = 1./(chi).^2 .* abs(runSum).^2;


Ea=Ee-Es;

end



function [stop]=endSum(A0, A1, num)
    stop=true;
    pDiff = abs((A1-A0)./A0) .* 100;
    
    for t=1:num,
        if((pDiff(t)>=0.001) || (A0(t) ==0))
            stop=false;
        end
    end

end
