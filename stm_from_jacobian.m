function [STM] = stm_from_jacobian(A, dt)
%STM_FROM_JACOBIAN
%
% Inputs:
%   A   - The Jacobian Matrix
%   dt  - Change in time

STM = expm(A*dt);
end