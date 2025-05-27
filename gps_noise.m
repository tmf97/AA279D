function sigma_pseudorange = gps_noise(r)
%GPS_NOISE 
%   Computes link budget parameters to give an estimate of the GPS
%   measurement noise given an r above earth.

c = 3e8; % m/s
bdll = 1/12; % Hz
gps_r = 26560000; % m
f_L1 = 1.57542e9;      % L1 frequency [Hz]

if r > gps_r
    distance_to_gps = gps_r + r;
    gps_antenna_gain = -10; % dBi
else
    distance_to_gps = gps_r - r;
    gps_antenna_gain = 3; % dBi
   
end
% link budget:

FSPL = 20*log10(distance_to_gps) + 20*log10(f_L1) + 20*log10(4*pi/c);
P_transmit = 48; % dBm

% RF GPS power received by sat
Gr = 16; % gain of the receiver on the spacecraft
Pr = P_transmit + gps_antenna_gain - FSPL + Gr;


% noise power
k = 1.38e-23;
T = 290;           % K
B = 2e6;           % 2 MHz, bandwidth of GPS L1 signal
N = 10*log10(k*T*B) + 30;

snr_db = Pr - N;
snr = 10.^(snr_db/10);
sigma_pseudorange = (c / 1.023e6) * sqrt(bdll./snr);
end