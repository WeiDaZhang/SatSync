close all;
satellite_start_delta = 49.65/180*pi;
satellite_start_alpha = 0;
satellite_altitude = 100e3; %(m)

ground_station_latitude = 50/180*pi;
ground_station_longitude = 0;
ground_station_height = 100;    % (m)

c = 299792458;      %m/s
average_earth_radius = 6.371e6;

% Ground Station Coordinates
[x0,y0,z0] = Geodetic2ECEF( ground_station_latitude, ...
                            ground_station_longitude, ...
                            ground_station_height);
figure(1);
plot3(x0,y0,z0,'o')
grid on;
hold on;

plot3(0,0,0,'+')

satellite_omega = SatelliteAngularVelocity(satellite_altitude + average_earth_radius);

% Signal Carrier frequency
f0 = 1.57542e9; %(Hz)
f0 = 1.584765e9; %(Hz)
lambda0 = c/f0;
w0 = 2*pi*f0;
k0 = w0/c;
p0dBm = 40; % (dBm)

% Oversampling Rate (How many samples each chip duration)
nOverSample = 12;
% Chip Rate (Rate of Spectrum Spread Code)
ChipRate = 4.095e6;
% Number of Point
N = 2^25;
% Local Time
t = 0 : 1/ChipRate/nOverSample : (N-1)/ChipRate/nOverSample;

% Generating Remote t'
PhaseNoiseVector = [1e-5 1e-5 1e-15 1e-10 1e-15];
t_prime = t + PhaseNoisePowerLaw2TIE(N, 1/ChipRate/nOverSample, PhaseNoiseVector);

% Satellite Coordinate at given time of sampling with starting orbital
% (angle) locations
[x,y,z] = ECS2ECEF( satellite_start_delta + t_prime*satellite_omega, ...
                    satellite_start_alpha, ...
                    satellite_altitude + average_earth_radius);

figure(1);
plot3(x,y,z)

% Satellite-Ground Station Range
sqr_R = ([(x - x0)',(y - y0)',(z - z0)']).^2;
range_x = sqrt(sum(sqr_R,2))';
% figure;plot(range_x)

% Received Power
link_loss = 20*log10(4*pi*range_x/lambda0);
received_power = p0dBm - link_loss;
received_amp = sqrt(10.^(received_power/10)*50);

% 1023 LFSR Polinomial x^10 + x^7 + 1
LFSR_POL = [12 11 10 4];
LFSR_INIT = [0 0 0 0 0 1];
ChipCode_b = LFSR_t(t_prime - range_x/c, ChipRate, LFSR_INIT, LFSR_POL);

received_r = received_amp.*sin(w0*t_prime - k0*range_x + ChipCode_b*pi);
figure;plot(abs(fft(received_r)))