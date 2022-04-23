% Klobuchar Ionospheric Model

function slant_factor = Klobuchar_Ionospheric_Model(latitude_user, longitude_user, elevation_sat, azimuth_sat )
% earth-centred angle
psi = 0.0137/(elevation_sat+0.11) - 0.022;

% latitude of the Ionospheric Pierce Point
latitude_IPP = latitude_user + psi * cos(azimuth_sat);
if(latitude_IPP > 0.416)
    latitude_IPP = 0.416;
elseif(latitude_IPP < -0.416)
    latitude_IPP = -0.416;
end

% longitude of the IPP
longitude_IPP = longitude_user + psi * sin(azimuth_sat) / cos(latitude_IPP);

% geomagnetic latitude of the IPP
latitude_m = latitude_IPP + 0.064 * cos(longitude_IPP - 1.617);

 % slant factor
 slant_factor = 1 + 16 * (0.53 - elevation_sat);