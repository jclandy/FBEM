function [ f , x , y ] = rsgene2D_anisotrop( N , M, rL , rW, h , Lx , Ly , LN)
%
% [f,x,y] = rsgene2D(N,rL,h,clx,cly) 
%
% generates a square 2-dimensional random rough surface f(x,y) with NxN 
% surface points. The surface has a Gaussian height distribution and 
% exponential autocovariance functions (in both x and y), where rL is the 
% length of the surface side, h is the RMS height and clx and cly are the 
% correlation lengths in x and y. 
%
% Input:    N   - number of surface points (along square side)
%           rL  - length of surface (along square side)
%           h   - rms height
%           L - correlation length (in x and y)
%
% Output:   f  - surface heights
%           x  - surface points
%           y  - surface points
%

% (C) Jack Landy & Alex Komarov, 2014 (adapted from code originally
% developed by David Bergström)

format long;

x = linspace(-rL/2,rL/2,N); y = linspace(-rW/2,rW/2,M);
[X,Y] = meshgrid(x,y);

if LN == 0
    h2 = h;
else
    h2 = 1;
end

Z = h2.*randn(M,N); % uncorrelated Gaussian random rough surface distribution
                                % with rms height h
                                
dx = rL / (N-1);
%ro = sqrt( X.^2+Y.^2 );

% Gaussian correlation function
%    C = exp( -ro.^2 / L ^2 );  % Gaussian correlation function
 %   F = 2*dx/(L*sqrt(pi)) * exp( - 2 * ro.^2 / L ^2 );   % precalculated filtering function for Gaussian correlation function

% Exponential correlation function
%L = Lx;

C = exp( - sqrt( (X / Lx).^2 + (Y / Ly).^2 ) );   % Exponential correlation function

%        C = exp( - ro / L );   % Exponential correlation function

%F =  dx / ( sqrt(2*pi) * L * gamma(0.75) ) * ( ( ro / (2*L) )  .^ ( -0.25 ) )  .* besselk( -0.25 , ro / L );                    % precalculated filtering function for Exponential correlation function

fft2_F = sqrt( fft2( C ) );

% correlation of surface including convolution (faltung), inverse
% Fourier transform and normalizing prefactors
fft2_Z = fft2(Z);
%f = 2*rL/N/clx*ifft2(fft2(Z).*fft2(F));
f = real(ifft2( fft2_Z.*fft2_F));     %using our method
%f2 = ifft2( fft2_Z.*fft2(F) );    % using analytically precalculated Filtering function

if LN == 1
    mu = log(1/sqrt(1 + h2^2/1^2));
    var = log(1 + h2^2/1^2);

    f = exp(mu + sqrt(var)*f);
    f = f*h*(1/std2(f));
    f = f - mean2(f);

else
end
    
end