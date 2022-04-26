% Geodetic Coordinates to Geocentric Cartesian Coordinates (ECEF)
% Geodetic coordinates are a type of curvilinear orthogonal coordinate system
% used in geodesy based on a reference ellipsoid.
% Geodetic different from Geocentric by considering earth ellipsoidal surface

function [x, y, z] = Geodetic2ECEF(latitude_phi, longitude_lambda, ellipsoidal_height)
% (Geodetic) latitude_phi (rad) is defined as the angle between the equatorial plane
% and the surface normal at a point on the ellipsoid.

% longitude_lambda (rad) measures the rotational angle between the zero meridian
% and the measured point.

% ellipsoidal_height (m), a.k.a. geodetic altitude is defined as the height
% above the ellipsoid surface, normal to the ellipsoid.



% check input length matching
vector_of_extend = convert_length([  length(latitude_phi) ...
                                     length(longitude_lambda) ...
                                     length(ellipsoidal_height)]);
if(vector_of_extend == 0)
    disp('length mishmatch.')
    return
end

% extend input to vector
if(vector_of_extend(1) ~= 1)
    latitude_phi = latitude_phi.*ones(1,vector_of_extend(1));
end
if(vector_of_extend(2) ~= 1)
    longitude_lambda = longitude_lambda.*ones(1,vector_of_extend(2));
end
if(vector_of_extend(3) ~= 1)
    ellipsoidal_height = ellipsoidal_height.*ones(1,vector_of_extend(3));
end


% earth_a (m) is the equatorial radius (semi-major axis).
earth_a = 6378137; % meter

% flattening factor earth_f = 1−earth_b/earth_a.
earth_f = 1/298.257223563;

% earth_b is the polar radius (semi-minor axis).
earth_b = (1 - earth_f)*earth_a;

% eccentricity sqr_earth_e = earth_e^2 is related with earth_a, earth_b and earth_f
sqr_earth_e = (earth_a^2 - earth_b^2)/earth_a^2;

% earth_N is the prime vertical radius of curvature, function of latitude_phi.
earth_N = earth_a./sqrt(1 - sqr_earth_e.*sin(latitude_phi).^2);

x = (earth_N + ellipsoidal_height).* cos(latitude_phi) .* cos(longitude_lambda);
y = (earth_N + ellipsoidal_height).* cos(latitude_phi) .* sin(longitude_lambda);
z = ((1 - sqr_earth_e) .* earth_N + ellipsoidal_height) .* sin(latitude_phi);
