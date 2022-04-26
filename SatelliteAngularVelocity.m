% Satellite Angular Velocity
function satellite_omega = SatelliteAngularVelocity(satellite_h)
% satellite_h (m) is the distance from the satellite to the centre earth
% satellite_omega (rad/s) Satellite orbital angular velocity

G = 6.673e-11; % Gravity Constant N*m^2/kg^2
M_earth = 5.98e24; % Earth Mass kg

satellite_velocity = sqrt(G * M_earth / satellite_h);

orbital_perimeter = 2 * pi * satellite_h;

satellite_omega = satellite_velocity/orbital_perimeter * 2 * pi;