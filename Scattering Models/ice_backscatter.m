function [theta,sigma_0_ice_surf,tau_ice] = ice_backscatter(lambda,sigma_si,l_si,T_si,S_si,h_s,beta_c,epsr_ds)

% Models the backscattering coefficient from the snow-sea ice interface,
% and interpolates the scattering signature to a spline function

% Input:
% lambda = radar wavelength, m
% sigma_si = small-scale rms height of interface, m
% l_si = small-scale correlation length of interface, m
% T_si = sea ice bulk temperature, C
% S_si = sea ice bulk salinity, psu
% h_s = snow depth, m (to calculate transmission coefficients)
% beta_c = Effective width of angular extent of coherent component
% epsr_ds = Relative permittivity of dry snow (if overlying snow layer
% present)

% Output:
% theta = angular sampling of the scattering signature
% sigma_0_ice_surf = smoothing spline characterizing the scattering
% signature of the snow-ice interface, dB
% tau_ice = transmission coefficient at snow-ice interface


% Uses the following codes from external sources:
% I2EM_Backscatter_model.m (Fawwaz Ulaby)
% RelDielConst_PureIce.m (Fawwaz Ulaby)
% TVBmodel_HeterogeneousMix.m (Fawwaz Ulaby)

% (C) Jack Landy, University of Bristol, 2018

%% Angular Sampling for Scattering Signature

theta = logspace(log10(1e-6),log10(pi/2),200);

%% Geophysical & Antenna Parameters

c = 299792458; % speed of light, m/s
f_c = c/lambda; % radar frequency, Hz
k0 = (2*pi)/lambda; % wavenumber

% mss_exp_si = (sqrt(2/pi)*sigma_si/l_si*sqrt(5*k*l_si - atan(5*k*l_si)))^2; % mean-square slope for exponential ACF (Dierking, 2000)

% Validity criteria for IEM
% max_sigma = 2/k; % maximum viable rms height for IEM
% 
% rayleigh = lambda/(8*cos(max_theta)); % surface = smooth if sigma < rayleigh criterion
% fraunhofer = lambda/(32*cos(max_theta)); % surface = smooth if sigma < fraunhofer criterion

% Pure ice density, kg/m^3
rho_i = 917 - 0.14*T_si; 

% Brine volume fraction & salinity from Cox and Weeks 1983, -2 to -23 C
F1_T = -4.732 + -2.245e1*T_si + -6.397e-1*T_si^2 + -1.074e-2*T_si^3;

V_b = (rho_i*1e-3*S_si)/F1_T; % brine volume
% T_f = -0.0575*S_si + 1.710523e-3*S_si^(3/2) - 2.154996e-4*S_si^2; % freezing T of seawater, C
S_br = 1000*(1 - 54.11/T_si)^-1; % brine salinity, in ppt

%% Dielectric Properties

% Sea ice dielectrics
[epsr_i,epsi_i] = RelDielConst_PureIce(T_si,f_c*1e-9); % permittvity of pure ice
[epsr_br,epsi_br] = RelDielConst_Brine(T_si,f_c*1e-9,S_br); % permittivity of brine inclusions

[eps_si] = TVBmodel_HeterogeneousMix(epsr_br + 1i*epsi_br , epsr_i + 1i*epsi_i, 2, V_b); % Mixture model with spherical inclusions

% Fresnel reflection & transmission coefficients
epsr_a = 1; % relative permittivity of air

if h_s > 0
    eta_0 = 376.73031346177; % intrinsic impedance of free space, ohms
    eta_1 = eta_0/sqrt(epsr_ds);
    eta_2 = eta_0/sqrt(real(eps_si));
    theta_2 = asin(sin(theta).*(sqrt(epsr_a/epsr_ds))); % transmission angle
else
    eta_0 = 376.73031346177; % intrinsic impedance of free space, ohms
    eta_1 = eta_0/sqrt(epsr_a);
    eta_2 = eta_0/sqrt(real(eps_si));
    theta_2 = asin(sin(theta).*(sqrt(epsr_a/real(eps_si)))); % transmission angle
end

rho_H = (eta_2*cos(theta) - eta_1*cos(theta_2))./(eta_2*cos(theta) + eta_1*cos(theta_2)); % reflection coeff
rho_V = (eta_2*cos(theta_2) - eta_1*cos(theta))./(eta_2*cos(theta_2) + eta_1*cos(theta));% reflection coeff

% rho_i_0 = (rho_H(1)+rho_V(1))/2; % ice reflection coeff at normal incidence

tau_H = 1 + rho_H; % transmission coeff
tau_V = (1 + rho_V)*(cos(theta)/cos(theta_2)); % transmission coeff
tau_ice = spline(theta,(tau_H + tau_V)/2);

gamma_H = rho_H.^2; % reflectivity (intensity)
gamma_V = rho_V.^2; % reflectivity (intensity)


%% Backscattering Coefficient of Snow-Ice Interface, sigma0

% Calculate coherent vs. incoherent surface scattering ratio
% psi = k0*sigma_si*cos(theta); % frequency-dependent roughness parameter
% omega = exp(-4*psi.^2); % fractional coherent component

% Calculate coherent reflected backscattering coefficient
% sigma_0_HH_coh = ((gamma_H.*omega)/beta_c^2)*exp(-4*k0^2*sigma_si^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient
% sigma_0_VV_coh = ((gamma_V.*omega)/beta_c^2)*exp(-4*k0^2*sigma_si^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient
sigma_0_HH_coh = 4*((gamma_H.*1)/beta_c^2)*exp(-4*k0^2*sigma_si^2).*exp(-4*theta.^2/beta_c^2); % equation 6 of Fung and Eom 1983
sigma_0_VV_coh = 4*((gamma_V.*1)/beta_c^2)*exp(-4*k0^2*sigma_si^2).*exp(-4*theta.^2/beta_c^2);

% Calculate incoherent surface backscattering coefficient
% Run single-scattering IEM for relevant range of facet incidence angles
sigma_0_HH_si_surf = zeros(length(theta),1);
sigma_0_VV_si_surf = zeros(length(theta),1);
for i = 1:length(theta)
    [sigma_0_VV_si_surf(i), sigma_0_HH_si_surf(i), ~] = I2EM_Backscatter_model(f_c*1e-9, sigma_si, l_si, theta(i)*180/pi, eps_si, 1, []);
end

% Calculate total co-polarized surface backscattering cofficients,
% including coherent reflected power
sigma_0_HH_si_surf = 10*log10(sigma_0_HH_coh + 10.^(sigma_0_HH_si_surf'/10)); % H-pol, dB
sigma_0_VV_si_surf = 10*log10(sigma_0_VV_coh + 10.^(sigma_0_VV_si_surf'/10)); % V-pol, dB

% Assuming no coherent reflected power
sigma_0_HH_si_surf(isinf(sigma_0_HH_si_surf))=NaN; sigma_0_VV_si_surf(isinf(sigma_0_VV_si_surf))=NaN;

% Build spline interpolants (assumption that scattering is polarization-independent)
% dB
sigma_0_ice_surf = spline(theta,(sigma_0_HH_si_surf + sigma_0_VV_si_surf)/2);


end

