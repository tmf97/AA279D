function [koes, times] = read_horizons(fname)
%INPUTS
%   fname   - (str) Full path to horizons output file
%    
%OUTPUTS
%   data     - Cell array containing each column from file
%
%NOTES
%
% Setup in Horizons is as follows:
%
% Ephemeris Type [change] : Osculating Orbital Elements
% Target Body [change] : 	MMS-#
% Coordinate Origin [change] : 	Earth (Geocentric) [500]
% Time Span [change] : 	Start=2015-10-17, Stop=2015-11-18, Step=1 h
% Table Settings [change] : CSV format=YES
% Display/Output [change] : download/save (plain text file)

if ~exist(fname,'file')
    disp([fname,' does not exist.']);
    return
end

fid = fopen(fname);
raw = textscan(fid, '%s', 'delimiter', '\n');
s = find(strcmpi(raw{1},'$$SOE'));
e = find(strcmpi(raw{1},'$$EOE'));
frewind(fid);

% Symbol meaning:
% 
%   JDTDB    Julian Day Number, Barycentric Dynamical Time
%     EC     Eccentricity, e
%     QR     Periapsis distance, q (km)
%     IN     Inclination w.r.t X-Y plane, i (degrees)
%     OM     Longitude of Ascending Node, OMEGA, (degrees)
%     W      Argument of Perifocus, w (degrees)
%     Tp     Time of periapsis (Julian Day Number)
%     N      Mean motion, n (degrees/sec)
%     MA     Mean anomaly, M (degrees)
%     TA     True anomaly, nu (degrees)
%     A      Semi-major axis, a (km)
%     AD     Apoapsis distance (km)
%     PR     Sidereal orbit period (sec)
% JDTDB Cal Date EC QR IN OM W  Tp N  MA TA A  AD PR
fstring = '%f %s %f %f %f %f %f %f %f %f %f %f %f %f';

data_raw = textscan(fid,fstring ,e-s+1,'delimiter', ',', 'headerlines', s);
fclose(fid);

% return only the stuff we care about, in a format that's useful, and with
% sensible units
KM_TO_M = 1000;

koes = [KM_TO_M * data_raw{12}, ... % semi-major axis
        data_raw{3}, ... % eccentricity
        deg2rad(data_raw{5}), ... % inclination
        deg2rad(data_raw{6}), ... % RAAN
        deg2rad(data_raw{7}), ... % argument of perigee
        deg2rad(data_raw{11})].'; % true anomaly

times = data_raw{1}.'; % in julian days

