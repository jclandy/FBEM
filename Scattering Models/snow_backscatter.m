function [theta,sigma_0_snow_surf,sigma_0_snow_vol,kappa_e,tau_snow,c_s,epsr_ds] = snow_backscatter(lambda,sigma_s,l_s,T_s,rho_s,r_s,h_s,beta_c)

% Models the backscattering coefficients from the snow-ice interface and
% snow volume, and interpolates scattering signatures to spline functions
% (Suitable only for dry snow)

% Input:
% lambda = radar wavelength, m
% sigma_s = small-scale rms height of interface, m
% l_s = small-scale correlation length of interface, m
% T_s = snow bulk temperature, C
% rho_s = snow bulk density, kg/m^3
% r_s = snow grain size, m
% h_s = snow depth, m
% beta_c = Effective width of angular extent of coherent component


% Output:
% theta = angular sampling of the scattering signature
% sigma_0_snow_surf = smoothing spline characterizing the scattering
% signature of the air-snow interface, dB
% sigma_0_snow_vol = smoothing spline characterizing the scattering
% signature of the snow volume, dB
% kappa_e = extinction coefficient of snow volume, Np/m
% tau_snow = transmission coefficient at air-snow interface
% c_s = speed of light in snow, m/s
% epsr_ds = relative permittivity of dry snow


% Uses the following codes from external sources:
% I2EM_Backscatter_model.m (Fawwaz Ulaby)
% RelDielConst_DrySnow.m (Fawwaz Ulaby)
% MieExtinc_DrySnow.m (Fawwaz Ulaby)

% (C) Jack Landy, University of Bristol, 2018

%% Angular Sampling for Scattering Signature

theta = logspace(log10(1e-6),log10(pi/2),200);

%% Geophysical & Antenna Parameters

c = 299792458; % speed of light, m/s
f_c = c/lambda; % radar frequency, Hz

% mss_exp_s = (sqrt(2/pi)*sigma_s/l_s*sqrt(5*k*l_s - atan(5*k*l_s)))^2; % mean-square slope for exponential ACF (Dierking, 2000)

% Reduced speed of light in snow, m/s
c_s = c*(1 + 0.51*rho_s*1e-3)^(-1.5); 


%% Dielectric Properties

% Dry snow dielectrics
[epsr_ds,epsi_ds] = RelDielConst_DrySnow(T_s,rho_s*1e-3,f_c*1e-9);
eps_ds = epsr_ds + 1i*epsi_ds;

% Fresnel reflection & transmission coefficients
epsr_a = 1; % relative permittivity of air

eta_0 = 376.73031346177; % intrinsic impedance of free space, ohms
eta_1 = eta_0/sqrt(epsr_a);
eta_2 = eta_0/sqrt(epsr_ds);
theta_2 = asin(sin(theta).*(sqrt(epsr_a/epsr_ds))); % transmission angle

rho_H = (eta_2*cos(theta) - eta_1*cos(theta_2))./(eta_2*cos(theta) + eta_1*cos(theta_2)); % reflection coeff
rho_V = (eta_2*cos(theta_2) - eta_1*cos(theta))./(eta_2*cos(theta_2) + eta_1*cos(theta));% reflection coeff

tau_H = 1 + rho_H; % transmission coeff, H-pol
tau_V = (1 + rho_V)*(cos(theta)/cos(theta_2)); % transmission coeff, V-pol
tau_snow = spline(theta,(tau_H + tau_V)/2);

% gamma_H = rho_H.^2; % reflectivity (intensity)
% gamma_V = rho_V.^2; % reflectivity (intensity)


%% Backscattering Coefficient of Air-Snow Interface, sigma0

% Calculate coherent vs. incoherent surface scattering ratio
% psi = k*sigma_s*cos(theta); % frequency-dependent roughness parameter
% omega = exp(-4*psi.^2); % fractional coherent component

% Calculate coherent reflected backscattering coefficient
% sigma_0_HH_coh = (gamma_H/beta_c^2)*exp(-4*k^2*sigma_s^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient
% sigma_0_VV_coh = (gamma_V/beta_c^2)*exp(-4*k^2*sigma_s^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient

% Calculate incoherent surface backscattering coefficient
% Run single-scattering IEM for relevant range of facet incidence angles
sigma_0_HH_s_surf = zeros(length(theta),1);
sigma_0_VV_s_surf = zeros(length(theta),1);
for i = 1:length(theta)
    [sigma_0_VV_s_surf(i), sigma_0_HH_s_surf(i), ~] = I2EM_Backscatter_model(f_c*1e-9, sigma_s, l_s, theta(i)*180/pi, eps_ds, 1, []);
end

% Calculate total co-polarized surface backscattering cofficients,
% including coherent reflected power
% sigma_0_HH_s_surf = 10*log10(sigma_0_HH_coh + 10.^(sigma_0_HH_s_surf'/10)); % H-pol, dB
% sigma_0_VV_s_surf = 10*log10(sigma_0_VV_coh + 10.^(sigma_0_VV_s_surf'/10)); % V-pol, dB

% Assuming no coherent reflected power
sigma_0_HH_s_surf(isinf(sigma_0_HH_s_surf))=NaN; sigma_0_VV_s_surf(isinf(sigma_0_VV_s_surf))=NaN;

% Build spline interpolants (assumption that scattering is polarization-independent)
% dB
sigma_0_snow_surf = spline(theta,(sigma_0_HH_s_surf + sigma_0_VV_s_surf)/2);


%% Backscattering Coefficient of Snow Volume, sigma0

% chi = sqrt(epsr_a)*(2*pi*r_s)/lambda; % normalized circumference
% n = sqrt(epsr_i/epsr_a);
% rayleigh_approximation = abs(n*chi); % Rayleigh scattering appropriate if <0.5

% Mie extinction coefficient in dry snow
[~,~,kappa_e,~] = MieExtinc_DrySnow(rho_s*1e-3,r_s,f_c*1e-9,T_s);

% Two-way loss factor
L_theta = exp((-2*kappa_e*h_s)./cos(theta_2)); % including scattering

% Penetration depth into snow cover (Ulaby et al 1984)
% Ignores scattering losses (valid when grain size <0.0001 m)
% delta_p = lambda/(4*pi)*(epsr_ds/2*(sqrt(1 + (epsi_ds/epsr_ds)^2) - 1))^(-1/2); % m
% 
% kappa_a2 = 1/delta_p; % Alternative absorption coeff calculated from penetration depth
% L_theta2 = exp((-2*kappa_a2*h_s)./cos(theta_2)); % ignoring scattering, following Nanden et al 2017

% Backscattered power from snow layer, following Winebrenner et al 1992
sigma_0_HH_s_vol = 10*log10(tau_H.^2.*(1 - L_theta));
sigma_0_VV_s_vol = 10*log10(tau_V.^2.*(1 - L_theta));

% Build spline interpolants (assumption that scattering is polarization-independent)
% dB
sigma_0_snow_vol = spline(theta,(sigma_0_HH_s_vol + sigma_0_VV_s_vol)/2);


end

