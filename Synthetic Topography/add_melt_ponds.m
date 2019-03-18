function [z_mp,surface_type_mp] = add_melt_ponds(z,surface_type,f_p)

%% Adds melt ponds with fraction f_p to surface topography

% Input:
% z = z coordinate of surface
% surface_type: 0 = lead/ocean, 1 = sea ice
% f_p = melt pond fraction

% Output:
% z_mp = new z coordinate of surface
% surface_type_mp: 0 = lead/ocean, 1 = sea ice, 2 = melt pond

% Cumulative distribution function
[f,z2] = ecdf(z(:));

% Identify height to fill f_p
mp_level = z2(find(f>=f_p,1));

% Fill to level
z_mp = z;
z_mp(z<=mp_level) = mp_level;

% Code ponds
surface_type_mp = surface_type;
surface_type_mp(z<=mp_level) = 2;

end

