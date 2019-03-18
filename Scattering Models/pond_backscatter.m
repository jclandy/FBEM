function [theta,sigma_0_mp_surf] = pond_backscatter(lambda,T_fw,beta_c,u_a)

% Models the backscattering coefficient from a melt pond surface, and
% interpolates the scattering signature to a spline function

% Input:
% lambda = radar wavelength, m
% T_fw = freshwater temperature, C
% beta_c = Effective width of angular extent of coherent component (TUNING PARAMETER)
% u_a = wind speed, m/s

% Output:
% theta = angular sampling of the scattering signature
% sigma_0_lead_surf = smoothing spline characterizing the scattering
% signature of the lead surface, dB


% Uses the following codes from external sources:
% RelDielConst_PureWater.m

% (C) Jack Landy, University of Bristol, 2018

warning('off','all')

%% Empirical relation between wind speed and pond roughness (Scharien et al. 2014)

sigma_mp = (0.982*u_a - 3.702)/1000;
l_mp = (2191*u_a^-2.621 + 12.79)/1000;

%% Angular Sampling for Scattering Signature

theta = logspace(log10(1e-6),log10(pi/2),200);

%% Antenna Parameters

c = 299792458; % speed of light, m/s
f_c = c/lambda; % radar frequency, Hz
k = (2*pi)/lambda; % wavenumber

% mss_exp_mp = (sqrt(2/pi)*sigma_mp/l_mp*sqrt(5*k*l_mp - atan(5*k*l_mp)))^2; % mean-square slope for exponential ACF (Dierking, 2000)

% Validity criteria for IEM
% max_sigma = 2/k; % maximum viable rms height for IEM
% 
% rayleigh = lambda/(8*cos(theta)); % surface = smooth if sigma < rayleigh criterion
% fraunhofer = lambda/(32*cos(theta)); % surface = smooth if sigma < fraunhofer criterion

%% Dielectric Properties

% Freshwater dielectrics
[epsr_fw,epsi_fw] = RelDielConst_PureWater(T_fw,f_c*1e-9); % permittivity of freshwater
eps_fw = epsr_fw + 1i*epsi_fw;

% Fresnel reflection coefficients
epsr_a = 1; % relative permittivity of air

eta_0 = 376.73031346177; % intrinsic impedance of free space, ohms
eta_1 = eta_0/sqrt(epsr_a);
eta_2 = eta_0/sqrt(real(eps_fw));
theta_2 = asin(sin(theta).*(sqrt(epsr_a/real(eps_fw)))); % transmission angle

rho_H = (eta_2*cos(theta) - eta_1*cos(theta_2))./(eta_2*cos(theta) + eta_1*cos(theta_2)); % reflection coeff
rho_V = (eta_2*cos(theta_2) - eta_1*cos(theta))./(eta_2*cos(theta_2) + eta_1*cos(theta));% reflection coeff

% tau_H = 1 + rho_H; % transmission coeff
% tau_V = (1 + rho_V)*(cos(theta)/cos(theta_2)); % transmission coeff

gamma_H = rho_H.^2; % reflectivity (intensity)
gamma_V = rho_V.^2; % reflectivity (intensity)


%% Backscattering Coefficient of Snow-Ice Interface, sigma0

% Calculate coherent vs. incoherent surface scattering ratio
psi = k*sigma_mp*cos(theta); % frequency-dependent roughness parameter
omega = exp(-4*psi.^2); % fractional coherent component

% Calculate coherent reflected backscattering coefficient (linear)
sigma_0_HH_coh = (gamma_H/beta_c^2)*exp(-4*k^2*sigma_mp^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient
sigma_0_VV_coh = (gamma_V/beta_c^2)*exp(-4*k^2*sigma_mp^2).*exp(-theta.^2/beta_c^2); % coherent component of backscattering coefficient

% Calculate incoherent surface backscattering coefficient
% Run single-scattering IEM for relevant range of facet incidence angles
sigma_0_HH_inc = zeros(length(theta),1);
sigma_0_VV_inc = zeros(length(theta),1);
for i = 1:length(theta)
    [sigma_0_VV_inc(i), sigma_0_HH_inc(i), ~] = I2EM_Backscatter_model(f_c*1e-9, sigma_mp, l_mp, theta(i)*180/pi, eps_fw, 1, []);
end

% Calculate total co-polarized surface backscattering cofficients,
% including coherent reflected power
sigma_0_HH_mp_surf = 10*log10(sigma_0_HH_coh + 10.^(sigma_0_HH_inc'/10).*(1-omega)); % H-pol, dB
sigma_0_VV_mp_surf = 10*log10(sigma_0_VV_coh + 10.^(sigma_0_VV_inc'/10).*(1-omega)); % V-pol, dB
sigma_0_HH_mp_surf(isinf(sigma_0_HH_mp_surf))=NaN; sigma_0_VV_mp_surf(isinf(sigma_0_VV_mp_surf))=NaN;

% Build spline interpolants (assumption that scattering is polarization-independent)
% dB
sigma_0_mp_surf = spline(theta,(sigma_0_HH_mp_surf + sigma_0_VV_mp_surf)/2);

warning('on','all')

end

