%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------------------------------------------------------- %%
%            SENSITIVITY DISTRIBUTION (REFLECTION COEFFICIENT)            %
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [SENS_IP,SENS_QP] = FDEM1DSENS_RC_cumulative(S,M,par)
%
%  Use:
%  Calculates the sensitivity distribution of a given layered soil medium 
%  and loop-loop configuration towards a certain physical property using 
%  the brute-force or perturbation method. Typical characteristics of the 
%  soil medium are stored in the Model structure (M) while the sensor 
%  characteristics are stored in the Sensor structure (S).
%
%  Input:
%  S (structure)           Sensor characteristics
%  M (structure)           Model characteristics
%  par                     Sensitivity parameter ('con','sus','perm')
%
%  Output:
%  SENS_IP                 IP sensitivity
%  SENS_QP                 QP sensitivity
%  Error (structure)       Estimated max. error (.IP and .QP)
%
%  [modifications]
%   JN, NOV 2019 - Remove the cumulative Sens as output
%   JN, OCT 2019 - SENS is output as ecdf
%   JN, OCT 2019 - original IP and QP are output arguments 
%
%  Created by Daan Hanssens
%  UGent, Belgium
%  January 2017
%
%  Cite:
%  Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M., 
%  and P. De Smedt. Frequency-Domain Electromagnetic Forward and 
%  Sensitivity Modeling: Practical Aspects of modeling a Magnetic Dipole 
%  in a Multilayered Half-Space. IEEE Geoscience and Remote Sensing 
%  Magazine, 7(1), 74-85
%

function [SENS_IP, SENS_QP, ErrorIP, ErrorQP,FWD_IP_ori, FWD_QP_ori]...
    = FDEM1DSENS_RC_cumulative(S,M,par)
%
% [SENS_IP, SENS_IP_cumu, SENS_QP, SENS_QP_cumu, ErrorIP, ErrorQP, FWD_IP_ori, FWD_QP_ori] 
% Was included as output the cumulative Sensitivity Analysis for IP and QP
%
    %% Allocate and initialize 
    FWD_IP_alt_p    = zeros(1, numel(M.(par)));
    FWD_QP_alt_p    = zeros(1, numel(M.(par)));
    FWD_IP_alt_n    = zeros(1, numel(M.(par)));
    FWD_QP_alt_n    = zeros(1, numel(M.(par)));
    
    %
    % Store original profile
    %
    
        op= M.(par);
        pert= op .* 0.01;    % Get relative perturbation (1%)
 
    %   
    % Get original response
    %

        [FWD_IP_ori,FWD_QP_ori] = FDEM1DFWD_RC(S,M);    
      
    %   
    % Calculate partial derivatives
    %
    
        % Loop Model layers
        for i= 1:numel(M.(par))

                
            %
            % Get altered response (forward)
            %
            
                M.(par)= op;                         % Get original profile
                M.(par)(i)= op(i) + pert(i);
                [FWD_IP_alt_p(i),FWD_QP_alt_p(i)] = FDEM1DFWD_RC(S,M);


            %
            % Get altered response (backward)
            %
            
                M.(par)(i)= op(i) - pert(i);
                [FWD_IP_alt_n(i),FWD_QP_alt_n(i)] = FDEM1DFWD_RC(S,M);

                
        end
        
    %
    % First derivative (Output)
    %
            
        SENS_QP= (FWD_QP_alt_p - FWD_QP_ori) ./ pert;
       	SENS_IP= (FWD_IP_alt_p - FWD_IP_ori) ./ pert;
          
                
    %
    % Second derivative
    %
            
    	SENS_QP_pert_sd= (FWD_QP_alt_p - 2.*FWD_QP_ori + ...
                                FWD_QP_alt_n) ./ pert .^2;
       	SENS_IP_pert_sd= (FWD_IP_alt_p - 2.*FWD_IP_ori + ...
                                FWD_IP_alt_n) ./ pert .^2;
    
    %
    % Estimate maximum error (Output)
    %
        ErrorQP = max(SENS_QP_pert_sd.* pert/2);
        ErrorIP = max(SENS_IP_pert_sd.* pert/2);

    %                             
    % Compute cumulative response
    % to be checked. last layer does not count for sensibility
    %                             
%         SENS_QP_cumu = cumsum(abs(SENS_QP(1,1:end-1)))./max(cumsum(abs(SENS_QP(1,1:end-1))));
%         SENS_IP_cumu = cumsum(abs(SENS_IP(1,1:end-1)))./max(cumsum(abs(SENS_IP(1,1:end-1))));
%         SENS_QP_cumu(1,end+1) = 1;
%         SENS_IP_cumu(1,end+1) = 1;
end
