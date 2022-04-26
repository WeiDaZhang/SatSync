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
x = 0;
y = 0;
z = 0;

% check input length matching
vector_of_extend = convert_length([length(delta) length(alpha) length(r)]);
if(vector_of_extend == 0)
    disp('length mishmatch.')
    return
end

% check distance input
if(r < r0)
    r = r + r0;
    disp('Input distance smaller than average earth radius.');
    disp('Input distance considered as max orbital height.');
end

% extend input to vector
if(vector_of_extend(1) ~= 1)
    delta = delta.*ones(1,vector_of_extend(1));
end
if(vector_of_extend(2) ~= 1)
    alpha = alpha.*ones(1,vector_of_extend(2));
end
if(vector_of_extend(3) ~= 1)
    r = r.*ones(1,vector_of_extend(3));
end

x = cos(delta).*cos(alpha).*r;
y = cos(delta).*sin(alpha).*r;
z = sin(delta).*r;
return