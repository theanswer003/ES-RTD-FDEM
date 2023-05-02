function [SENS_IP_sus, SENS_QP_sus, SENS_IP_con, SENS_QP_con, FWD_IP, FWD_QP,...
    ErrorConIp, ErrorConQp, ErrorSusIp, ErrorSusQp] =...
    FDEM1DSENS_RC_2SENS_array(i,j,ECmodel,MSmodel,S,M)
%% --------------------------------------------------------------------- %%
%          BUILD ONE ARRAY WITH THE FORWARD RESPONSE (IP AND QP)          %
%             AND THE SENSITIVITY ANALYSIS OF BOTH EC AND MS              %
%               USING THE REFLECTION COEFFICIENT APPROACH                 %
% ----------------------------------------------------------------------- %
%
% Input:
%
% S - sensor characteristics
% M - model characteristics
% ECmodel - Previously simulated model of Electrical Conductivity.
% MSmodel - Previously simulated model of Magnetic Susceptibility.
%
% Output:
%
% FWD_IP - In-Phase (or Real) response (in ppm) of the magnetic field 
% FWD_QP - Quadrature-Phase (or Imaginary) response (in ppm) of the
%          magnetic field.
% SENS_IP_sus - In-Phase (or Real) sensitivity analysis using MS.
% SENS_QP_sus - Quadrature-Phase (or Imaginary) analysis analysis using MS.
% SENS_IP_con - In-Phase (or Real) sensitivity analysis using EC. 
% SENS_QP_con - Quadrature-Phase (or Imaginary) analysis analysis using EC.
%
%
% Functions:
%
% FDEM1DSENS_RC_cumulative - Calculates the sensitivity distribution of a 
%                            given layered soil medium and loop-loop
%                            configuration towards a certain physical 
%                            property using the brute-force or perturbation 
%                            method.
%
% MODIFICATIONS:
%
%
% João Narciso and Leonardo Azevedo
% CERENA/IST, Portugal
% October 2019
%

M.con = squeeze(ECmodel(i ,j , :))';  

M.sus = squeeze(MSmodel(i ,j , :))';              

[SENS_IP_sus, SENS_QP_sus, ErrorSusIp, ErrorSusQp, FWD_IP, FWD_QP] = FDEM1DSENS_RC_cumulative(S,M,'sus');   % Reflection Coefficient approach 
[SENS_IP_con, SENS_QP_con, ErrorConIp, ErrorConQp]                 = FDEM1DSENS_RC_cumulative(S,M,'con');   % Reflection Coefficient approach 

end