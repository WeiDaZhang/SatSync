% Equatorial Coordinate System to ECEF
% The positions of artificial Earth satellites are specified in geocentric 
% equatorial coordinates, also known as geocentric equatorial inertial (GEI),
% Earth-centred inertial (ECI), and conventional inertial system (CIS),
% all of which are equivalent in definition to
% the astronomical geocentric equatorial rectangular frames, above.
function [x, y, z] = ECS2ECEF(delta, alpha, r)
% delta (rad) is the angular distance of an object perpendicular to the celestial equator,
% positive to the north, negative to the south.

% alpha (rad) is  angular distance of an object eastward along the celestial equator
% from the vernal equinox to the hour circle (great circle) passing through the object.

% r (m) is the distance from the object to the centre earth

r0 = 6.371e6;   % average earth radius in meter
if(r < r0)
    r = r + r0;
    disp('Input distance smaller than average earth radius.');
    disp('Input distance considered as max orbital height.');
end

x = cos(delta)*cos(alpha)*r;
y = cos(delta)*sin(alpha)*r;
z = sin(delta)*r;
return