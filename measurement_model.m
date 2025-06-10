function [state_measurement, chief_measurement, deputy_measurement]  = measurement_model(pv_chief,pv_deputy, sigma)
%MEASUREMENT_MODEL GPS Measurement model
%   
% Inputs:
%   state_estimate  - state vector
%   sigma           - standard deviation of noise
%
% Output:
%   measurement_estimate     - estimate of measurement


noise_chief = sigma.* randn(6,1) .* [1, 1, 1, 0.1, 0.1, 0.1].';
noise_deputy = sigma.* randn(6,1) .* [1, 1, 1, 0.1, 0.1, 0.1].';
chief_measurement = pv_chief + noise_chief;
deputy_measurement = pv_deputy + noise_deputy;
state_measurement = pv2roe(chief_measurement, deputy_measurement);
end