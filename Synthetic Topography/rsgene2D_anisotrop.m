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
    mu_norm = 0;
    sigma_norm = h;
else
    % calculate required normal parameters, mean = 1
    mu_norm = log(1^2/sqrt(h^2 + 1^2));
    sigma_norm = sqrt(log(h^2/1^2 + 1));
end

Z = sigma_norm.*randn(M,N) + mu_norm; % uncorrelated Gaussian random rough surface distribution                         % with rms height h
                                
dx = rL / (N-1);

% Exponential correlation function
C = exp( - sqrt( (X / Lx).^2 + (Y / Ly).^2 ) );  
fft2_F = sqrt( fft2( C ) );

% correlation of surface including convolution (faltung), inverse
% Fourier transform and normalizing prefactors
fft2_Z = fft2(Z);
f = real(ifft2( fft2_Z.*fft2_F));     %using our method

% scale to original parameter values
f = f*sigma_norm/std(f(:));
f = f - (mean(f(:)) - mu_norm);

if LN == 1
    f = exp(f) - 1;
else
end
    
end