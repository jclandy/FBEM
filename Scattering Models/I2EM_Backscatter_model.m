%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.1: I2EM Backscattering from Single-Scale Random Surface

%Description: Code computes sigma_0_vv, sigma_0_hh, and
%sigma_0_hv for single-scale random surface with a specified correlation
%function

%Input Variables: 
    %er: complex dielectric constant of the scattering medium
    %thi: Incidence angle (deg)
    
    %sp: type of correlation function: 1- exponential, 2- Gaussian, 3-
            %x-power
    %xx: coefficient for the x-power correlation function 
    %sig: rms height (m)
    %L: correlation length (m)
    %fr: frequency (GHz)
    
%Output Products:
    % sigma_0_vv, sigma_0_hh, and sigma_0_hv in dB

%Book Reference: Section 10-3.9

%Matlab Code: 

function [sigma_0_vv, sigma_0_hh, sigma_0_hv] = I2EM_Backscatter_model(...
    fr, sig, L, thi, er, sp, xx)


%--the code calls two other functions. The first is the I2EM_Bistatic_model
%code which is used to calculate the co-polarized responses in backscatter.
%The second is the IEMX_model used to calculate the cross-polarized
%response. The later code is slow as it involved double integrations.


%--- The co-pol components:
ths = thi;
phs = 180;
[sigma_0_vv, sigma_0_hh] = I2EM_Bistat_model(fr, sig, L, thi, ths, phs, er, sp, xx);


%--- The cross-pol component
auto = 1; % auto = 1 allows for automatic selection of the number of spectral components 
          % auto = 0 forces the number of spectral components to 15 always.
          % selection of auto = 1 results in a slower code

[sigma_0_hv] = IEMX_model(fr, sig, L, thi, er, sp,xx, auto);
 
end


%Code 10.3: I2EM Bistatic Scattering from Single-Scale Random Surface

%Description: Code computes sigma_0_vv(thi, ths,phs) and
%sigma_0_hh(thi, ths,phs) for single-scale random surface

%Input Variables: 
    %er: complex dielectric constant of the scattering medium
    %thi: Incidence angle (deg)
    %ths: Scattering angle (deg)
    %phs: Relative azimuth angle (deg)
    
    %sp: type of correlation function: 1- exponential, 2- Gaussian
    %       3- x-power,  4- x-exponential
    %xx: coefficient (>1) needed for x-power and x-exponential correl. fnc.
    %sig: rms height (m)
    %L: correlation length (m)
    %fr: frequency (GHz)
    
%Output Products:
    % sigma_0_vv, sigma_0_hh in dB

%Book Reference: Section 10-3.9

%Matlab Code: 

function [sigma_0_vv, sigma_0_hh] = I2EM_Bistat_model(fr, sig, L, thi, ths, phs, er, sp, xx)


error = 1.0e8;
sig = sig * 100; % change from m to cm scale
L = L * 100;
mu_r = 1; % relative permeability

k = 2*pi *fr/30; % wavenumber in free space. Speed of light is in cm/sec
theta = thi .*pi/180; % transform to radian
phi = 0;
thetas = ths * pi/180;
phis = phs * pi/180;

ks = k * sig; % roughness parameter
kl = k * L;  

ks2 = ks .* ks; 

cs = cos(theta+ 0.01);
s = sin(theta+ 0.01);

sf = sin(phi);
cf = cos(phi);

ss = sin(thetas);
css = cos(thetas);

cfs = cos(phis);
sfs = sin(phis);

s2 = s * s;
% sq = sqrt(er - s2);

kx = k .* s .*cf;
ky = k .* s .*sf;
kz = k .* cs;

ksx = k .* ss .*cfs;
ksy = k .* ss .*sfs;
ksz = k .* css;

%-- reflection coefficients
rt = sqrt(er - s2);
Rvi = (er .*cs - rt) ./(er.*cs +rt);
Rhi = (cs - rt)./(cs + rt);


wvnb = k .* sqrt( (ss .*cfs - s .*cf).^2 + (ss .* sfs - s .* sf).^2 );

Ts = 1;

while error > 1.0e-8,
    Ts = Ts + 1;
    error = (ks2 .*(cs + css).^2 ).^Ts ./ factorial(Ts); 
end

%---------------- calculating roughness spectrum -----------

[wn, rss] = roughness_spectrum(sp,xx, wvnb, sig, L, Ts);


%----------- compute R- transition ------------

Rv0 = (sqrt(er)-1) ./(sqrt(er)+1);
Rh0 = -Rv0;

Ft = 8 * Rv0.^2 * ss *(cs + sqrt(er - s2))./(cs .* sqrt(er - s2));
a1 = 0;
b1 = 0;
for n = 1:Ts
    a0 = (ks .*cs).^(2*n) ./factorial(n);
    a1 = a1 + a0 *wn(n);
    b1 = b1 + a0 * (abs(Ft./2 + 2.^(n+1) .*Rv0./cs .*exp(-(ks.*cs).^2))).^2 ...
        * wn(n);
end
St = 0.25 * (abs(Ft)).^2 * a1 ./ b1;

St0 = 1 ./ (abs(1 + 8 *Rv0./(cs .* Ft))).^2;

Tf = 1 - St ./St0;


%----------- compute average reflection coefficients ------------
%-- these coefficients account for slope effects, especially near the
%brewster angle. They are not important if the slope is small.

sigx = 1.1 .*sig/L;
sigy = sigx;
xxx = 3*sigx;

Rav = dblquad(@(Zx, Zy)Rav_integration(Zx, Zy, cs,s,er,s2,sigx, sigy),-xxx,xxx, -xxx, xxx );

Rah = dblquad(@(Zx, Zy)Rah_integration(Zx, Zy, cs,s,er,s2,sigx, sigy),-xxx,xxx, -xxx, xxx );

Rav = Rav ./(2*pi * sigx * sigy);
Rah = Rah ./(2*pi * sigx * sigy);



%-- select proper reflection coefficients

if thi == ths && phs==180, %i.e. operating in backscatter mode
    Rvt = Rvi + (Rv0 - Rvi) .*Tf;
    Rht = Rhi + (Rh0 - Rhi) .*Tf;

else        % in this case, it is the bistatic configuration and average R is used
%     Rvt = Rav + (Rv0 - Rav) .* Tf;
%     Rht = Rah + (Rh0 - Rah) .* Tf;
    Rvt = Rav;
    Rht = Rah;
end

fvv = 2 .* Rvt .*(s .* ss - (1 + cs .* css) .* cfs)./(cs + css);
fhh = -2 .* Rht .*(s .* ss - (1 + cs .* css) .* cfs)./(cs + css);

%------- Calculate the Fppup(dn) i(s) coefficients ----
[Fvvupi, Fhhupi] = Fppupdn_is_calculations(+1, 1, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs);
[Fvvups, Fhhups] = Fppupdn_is_calculations(+1, 2, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs);
[Fvvdni, Fhhdni] = Fppupdn_is_calculations(-1, 1, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs);
[Fvvdns, Fhhdns] = Fppupdn_is_calculations(-1, 2, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs);


qi = k .* cs;
qs = k .* css;

%----- calculating  Ivv and Ihh ----

Ivv = zeros(Ts, 1); Ihh = Ivv;
    
for n = 1:Ts
    Ivv(n) = (kz + ksz).^n .* fvv .* exp(-sig^2 .* kz .* ksz) + ...
        0.25*(Fvvupi .*(ksz-qi).^(n-1) .*exp(-sig^2 .*(qi.^2 - qi.*(ksz-kz)))+ ...
        Fvvdni .*(ksz+qi).^(n-1) .*exp(-sig^2 .*(qi.^2 + qi.*(ksz-kz)))+ ...
        Fvvups .*(kz+qs).^(n-1) .*exp(-sig^2 .*(qs.^2 - qs.*(ksz-kz)))+ ...
        Fvvdns .*(kz-qs).^(n-1) .*exp(-sig^2 .*(qs.^2 + qs.*(ksz-kz))));

    
    Ihh(n) = (kz + ksz).^n .* fhh .* exp(-sig^2 .* kz .* ksz) + ...
        0.25*(Fhhupi .*(ksz-qi).^(n-1) .*exp(-sig^2 .*(qi.^2 - qi.*(ksz-kz)))+ ...
        Fhhdni .*(ksz+qi).^(n-1) .*exp(-sig^2 .*(qi.^2 + qi.*(ksz-kz)))+ ...
        Fhhups .*(kz+qs).^(n-1) .*exp(-sig^2 .*(qs.^2 - qs.*(ksz-kz)))+ ...
        Fhhdns .*(kz-qs).^(n-1) .*exp(-sig^2 .*(qs.^2 + qs.*(ksz-kz))));
end


%-- Shadowing function calculations

if thi==ths && phs==180 %i.e. working in backscatter mode
    ct = cot(theta);
    cts = cot(thetas);
    rslp = rss;
    ctorslp = ct / sqrt(2) ./rslp;
    ctsorslp = cts / sqrt(2) ./rslp;
    shadf = 0.5 *(exp(-ctorslp.^2) ./ sqrt(pi)./ctorslp - erfc(ctorslp));
    shadfs = 0.5 *(exp(-ctsorslp.^2) ./ sqrt(pi)./ctsorslp - erfc(ctsorslp));
    ShdwS = 1./(1 + shadf + shadfs); 
else
 ShdwS = 1;
end

%------- calculate the values of sigma_note --------------

sigmavv = 0;sigmahh =0; 

for n = 1:Ts
    a0 = wn(n) ./factorial(n) .*sig.^(2*n);
    
    sigmavv = sigmavv+ abs(Ivv(n)).^2 .*a0;
    sigmahh = sigmahh+ abs(Ihh(n)).^2 .*a0;

end

sigmavv = sigmavv * ShdwS * k^2 ./2 * exp(-sig.^2 .*(kz.^2 +ksz.^2));  
sigmahh = sigmahh * ShdwS * k^2 ./2 * exp(-sig.^2 .*(kz.^2 +ksz.^2));  


ssv = 10 * log10(sigmavv); 
ssh = 10 * log10(sigmahh); 

sigma_0_vv = ssv;
sigma_0_hh = ssh;

end

function [vv, hh] = Fppupdn_is_calculations(ud, is, Rvi,Rhi,er,k,kz,ksz,s,cs,ss,css,cf,cfs,sfs)

if is==1
    Gqi = ud .* kz;
    Gqti = ud .*k .*sqrt(er-s.^2);
    qi = ud .* kz;
    
    c11 = k .* cfs .*(ksz - qi);
    c21 = cs .*(cfs .*(k^2 .*s.*cf.*(ss .*cfs - s .* cf) + Gqi.*(k .* css - qi)) ...
        + k^2 .*cf .* s .*ss .*sfs^2);
    c31 = k.*s.*(s.*cf.*cfs.*(k.*css-qi) - Gqi.*(cfs.*(ss.*cfs -s.*cf)+ ss .*sfs^2));
    c41 = k .*cs.*(cfs.*css.*(k.*css - qi) + k .*ss.*(ss.*cfs-s.*cf));
    c51 = Gqi.*(cfs .*css.*(qi-k.*css) - k .*ss.*(ss.*cfs-s.*cf));
    
    c12 = k .* cfs .*(ksz - qi);
    c22 = cs .*(cfs .*(k^2 .*s.*cf.*(ss .*cfs - s .* cf) + Gqti.*(k .* css - qi)) ...
        + k^2 .*cf .* s .*ss .*sfs^2);
    c32 = k.*s.*(s.*cf.*cfs.*(k.*css-qi) - Gqti.*(cfs.*(ss.*cfs -s.*cf)- ss .*sfs^2));
    c42 = k .*cs.*(cfs.*css.*(k.*css - qi) + k .*ss.*(ss.*cfs-s.*cf));
    c52 = Gqti.*(cfs .*css.*(qi-k.*css) - k .*ss.*(ss.*cfs-s.*cf));    
end

if is==2
    Gqs = ud .* ksz;
    Gqts = ud .*k .*sqrt(er-ss.^2);
    qs = ud .* ksz;
    
    c11 = k .* cfs .*(kz + qs);
    c21 = Gqs .*(cfs.*(cs.*(k.*cs+qs)-k.*s.*(ss .*cfs-s.*cf))-k.*s.*ss.*sfs.^2);
    c31 = k .*ss.*(k.*cs.*(ss.*cfs - s.*cf)+ s.*(kz+qs));
    c41 = k.*css.*(cfs.*(cs.*(kz+qs)-k.*s.*(ss.*cfs-s.*cf))-k.*s.*ss.*sfs.^2);
    c51 = -css .*(k.^2 .*ss .*(ss.*cfs -s.*cf)+ Gqs.*cfs.*(kz+qs));
         
    c12 = k .* cfs .*(kz + qs);
    c22 = Gqts .*(cfs.*(cs.*(kz+qs)-k.*s.*(ss .*cfs-s.*cf))-k.*s.*ss.*sfs.^2);
    c32 = k .*ss.*(k.*cs.*(ss.*cfs - s.*cf)+ s.*(kz+qs));
    c42 = k.*css.*(cfs.*(cs.*(kz+qs)-k.*s.*(ss.*cfs-s.*cf))-k.*s.*ss.*sfs.^2);
    c52 = -css .*(k.^2 .*ss .*(ss.*cfs -s.*cf)+ Gqts.*cfs.*(kz+qs));
end

q=kz;
qt = k .*sqrt(er - s.^2);

vv = (1+Rvi) .*( -(1-Rvi) .*c11 ./q + (1+Rvi) .*c12 ./ qt) + ...
    (1 - Rvi) .*( (1-Rvi) .*c21 ./q - (1+Rvi) .*c22 ./ qt) + ...
    (1+Rvi) .*( (1-Rvi) .*c31 ./q - (1+Rvi) .*c32 ./er ./qt) + ...
    (1 - Rvi) .*( (1+Rvi) .*c41 ./q - er.*(1 - Rvi) .*c42 ./ qt) + ...
    (1+Rvi) .*( (1+Rvi) .*c51 ./q - (1-Rvi) .*c52 ./ qt);

hh = (1 + Rhi) .*( (1-Rhi) .*c11 ./q - er.*(1+Rhi) .*c12 ./ qt) - ...
    (1 - Rhi) .*( (1-Rhi) .*c21 ./q - (1+Rhi) .*c22 ./ qt) - ...
    (1 + Rhi) .*( (1-Rhi) .*c31 ./q - (1+Rhi) .*c32 ./qt) - ...
    (1 - Rhi) .*( (1+Rhi) .*c41 ./q - (1 - Rhi) .*c42 ./ qt) - ...
    (1 + Rhi) .*( (1+Rhi) .*c51 ./q - (1-Rhi) .*c52 ./ qt);

end


function Rav = Rav_integration(Zx, Zy, cs,s,er,s2,sigx, sigy)

A = cs + Zx .* s;
B = er .* (1 + Zx.^2 + Zy.^2);
CC = s2 - 2.*Zx .*s .*cs + Zx.^2 .* cs^2 + Zy.^2;

Rv = (er.*A - sqrt(B-CC))./(er.*A + sqrt(B-CC)); 

pd = exp(-Zx.^2 ./(2*sigx.^2) -Zy.^2 ./(2*sigy.^2));
Rav = Rv .* pd;

end

function Rah = Rah_integration(Zx, Zy, cs,s,er,s2,sigx, sigy)

A = cs + Zx * s;
B = er .* (1 + Zx.^2 + Zy.^2);
CC = s2 - 2.*Zx .*s .*cs + Zx.^2 .* cs^2 + Zy.^2;

Rh = (A - sqrt(B-CC))./(A + sqrt(B-CC)); 

pd = exp(-Zx.^2./(2*sigx.^2) -Zy.^2./(2*sigy.^2));
Rah = Rh .* pd;

end

function [wn, rss] = roughness_spectrum(sp,xx, wvnb, sig, L, Ts)

wn = zeros(Ts,1);
%-- exponential correl func
if sp ==1,
    for n = 1: Ts
        wn(n) = L.^2 / n.^2 .* (1+(wvnb.*L/n).^2).^(-1.5);
    end
    rss= sig ./L;
end

%-- gaussian correl func
if sp ==2,
    for n = 1: Ts
        wn(n) =  L.^2/(2*n) .* exp(-(wvnb.*L).^2 ./(4*n));
    end
    rss = sqrt(2) * sig./L;
end

%-- x-power correl func
if sp == 3,
   
    for n =1: Ts
      if wvnb==0, 
          wn(n) = L.^2 ./(3 * n - 2);
      else
          wn(n) = L.^2 * (wvnb.*L).^(-1 + xx*n) .* besselk(1 -xx*n, wvnb*L) ...
              ./ (2.^(xx*n -1) .* gamma(xx .*n)); 
      end
    end
%    if xx == 1.5
        rss  = sqrt(xx *2) .*sig ./L;
%    else
%        rss = 0;
%    end
end

%-- x- exponential correl func

if sp == 4,
    for n =1: Ts
        tmp = quad(@(z) x_exponential_spectrum(z,wvnb,L,n,xx), 0,9);
        wn(n) = L.^2 ./ n.^(2/xx) .* tmp;
    end
    rss  =  sig ./ L;
end

end

function [tmp] = x_exponential_spectrum(z,wvnb,L,n,xx)

tmp = exp(-abs(z).^xx) .* besselj(0, z*wvnb*L./(n.^(1/xx))).*z;
end


function sigvh = IEMX_model(fr, sig, L, theta_d, er, sp,xx, auto)

sig = sig * 100 ; % change to cm scale
L = L * 100; % change to cm scale;


%- fr: frequency in GHz
%- sig: rms height of surface in cm
%- L: correlation length of surface in cm
%- theta_d: incidence angle in degrees
%- er: relative permittivity
%- sp: type of surface correlation function

error = 1.0e8;

k = 2*pi *fr/30; % wavenumber in free space. Speed of light is in cm/sec
theta = theta_d .*pi/180; % transform to radian

ks = k * sig; % roughness parameter
kl = k * L;  

ks2 = ks .* ks; 
kl2 = kl.^2;

cs = cos(theta);
s = sin(theta+ 0.001);

s2 = s.^2;

%-- calculation of reflection coefficints
rt = sqrt(er - s2);

rv = (er *cs - rt) ./(er*cs +rt);
rh = (cs - rt)./(cs + rt);

rvh = (rv - rh) ./2;

%-- rms slope values
sig_l = sig/L;
if sp==1                %-- exponential correl func
    rss = sig_l;
end
if sp==2                %-- Gaussian correl func
    rss = sig_l * sqrt(2);
end

if sp==3                 %-- 1.5-power spectra correl func
    rss = sig_l * sqrt(2*xx);
end

%--- Selecting number of spectral components of the surface roughness
if auto == 0 
n_spec = 15; % number of terms to include in the surface roughness spectra
end

if auto == 1
   n_spec = 1; 
   while error > 1.0e-8,
    n_spec = n_spec + 1;
    error = (ks2 .*( 2*cs).^2 ).^n_spec ./ factorial(n_spec); 
   end

end

%-- calculating shadow consideration in single scat (Smith, 1967)

ct = cot(theta+ 0.001);
farg = ct /sqrt(2) ./rss;
gamma = 0.5 *(exp(-farg.^2) / 1.772 / farg - erfc(farg));
Shdw = 1 ./ (1 + gamma);

%-- calculating multiple scattering contribution
%------ a double integration function

svh = dblquad(@(r,phi)xpol_integralfunc(r, phi, sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec), 0.1, 1, 0, pi);

sigvh = 10*log10(svh .* Shdw);

end

function y = xpol_integralfunc(r, phi, sp, xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec)

cs2 = cs .^2;

r2 = r.^2;
nr = length(r);

sf = sin(phi);
csf = cos(phi);
rx = r .* csf;
ry = r .* sf;

%-- calculation of the field coefficients
rp = 1 + rvh;
rm = 1 - rvh; 

q = sqrt(1.0001 - r2);
qt = sqrt(er - r2);

a = rp ./q;
b = rm ./q;
c = rp ./qt;
d = rm ./qt;

%--calculate cross-pol coefficient
B3 = rx .* ry ./cs;
fvh1 = (b-c).*(1- 3*rvh) - (b - c./er) .* rp; 
fvh2 = (a-d).*(1+ 3*rvh) - (a - d.*er) .* rm;
Fvh = ( abs( (fvh1 + fvh2) .*B3)).^2;


%-- calculate shadowing func for multiple scattering 
au = q ./r ./1.414 ./rss;
fsh = (0.2821./au) .*exp(-au.^2) -0.5 .*(1- erf(au));
sha = 1./(1 + fsh); 

%-- calculate expressions for the surface spectra
wn = spectrm1(sp, xx, kl2, L, rx, ry, s, n_spec, nr);
wm = spectrm2(sp, xx, kl2, L, rx, ry, s, n_spec, nr);


%--compute VH scattering coefficient
acc = exp(-2* ks2 .*cs2) ./(16 .* pi);
vhmnsum = zeros(1,nr);
for n = 1: n_spec
    for m = 1: n_spec
        vhmnsum = vhmnsum + wn(n,:).*wm(m,:) .*(ks2*cs2).^(n+m) ...
            ./factorial(n)./factorial(m); 
    end
end

VH = 4 * acc .* Fvh .* vhmnsum .*r;
y = VH .* sha;

end

function wn = spectrm1(sp, xx, kl2, L, rx, ry, s, np, nr)

wn = zeros(np,nr);

if sp == 1  % exponential
    for n = 1: np       
        wn(n, :) = n* kl2 ./(n.^2 + kl2 *((rx-s).^2+ry.^2)).^1.5;
    end
end

if sp == 2  %  gaussian
    for n = 1: np
        wn(n,:) = 0.5 * kl2 ./n .* exp(-kl2*((rx-s).^2 + ry.^2)/(4*n)) ;
    end
end

if sp== 3 % x-power
    for n = 1: np
        wn(n,:) = kl2 ./(2.^(xx*n-1) .*gamma(xx*n)).* ( ( (rx-s).^2 ...
          + ry.^2).*L).^(xx*n-1) .* besselk(-xx*n+1, L*((rx-s).^2 + ry.^2));
    end
end

end

function wm = spectrm2(sp, xx, kl2, L, rx, ry, s, np, nr)

wm = zeros(np,nr);

if sp == 1  % exponential
    for n = 1: np
        wm(n,:) = n* kl2 ./(n.^2 + kl2 *((rx+s).^2+ry.^2)).^1.5;
    end
end

if sp == 2  %  gaussian
    for n = 1: np
        wm(n,:) = 0.5 * kl2 ./n .* exp(-kl2*((rx+s).^2 + ry.^2)/(4*n)) ;
    end
end

if sp== 3 % x-power
    for n = 1: np
        wm(n,:) = kl2 ./(2.^(xx*n-1) .*gamma(xx*n)).* ( ( (rx+s).^2 ...
          + ry.^2).*L).^(xx*n-1) .* besselk(-xx*n+1, L*((rx+s).^2 + ry.^2));
    end
end
end


