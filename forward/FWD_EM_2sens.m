function [FWD_IP, FWD_QP, SENS_IP_Sus, SENS_QP_Sus, SENS_IP_Con, SENS_QP_Con,...
    ErrorConIp, ErrorConQp, ErrorSusIp, ErrorSusQp] = FWD_EM_2sens(S,M, ECmodel, MSmodel, par, EMapprox)
%% --------------------------------------------------------------------- %%
%   FORWARD RESPONSE AND SENSITIVITY ANALYSIS OF THE SIMULATED EM MODEL   %
% ----------------------------------------------------------------------- %
%
% Computes the forward response (ppm) and the sensitivity distribution of
% a given layered half-space and loop-loop configuration. 
% 
%
% Input:
%
% S - sensor characteristics
% M - model characteristics
% ECmodel - Previously simulated model of Electrical Conductivity.
% MSmodel - Previously simulated model of Magnetic Susceptibility.
% par - parameter used for the calculation of the sensitivity of the model
%       (electrical conductivity / magnetic susceptibility / dielectric
%       permitivity)
% EMapprox - type of approach in the forward model 
%            (PM = Propagation Matrix approach;
%            (RM = Reflection Coefficient approach)
%
% Output:
%
% FWD_IP - In-Phase (or Real) response (in ppm) of the magnetic field 
% FWD_QP - Quadrature-Phase (or Imaginary) response (in ppm) of the
%          magnetic field.
% SENS_IP_sus - Sensitivity analysis of In-Phase (or Real) toward changes in MS.
% SENS_QP_sus - Sensitivity analysis of Quadrature-Phase (or Imaginary) toward changes in MS.
% SENS_IP_con - Sensitivity analysis of In-Phase (or Real) toward changes in EC. 
% SENS_QP_con - Sensitivity analysis of Quadrature-Phase (or Imaginary) toward changes in EC.
% SENS_IP_sus_cumu - Cumulative sensitivity analysis of In-Phase (or Real) toward changes in MS.
% SENS_QP_sus_cumu - Cumulative sensitivity analysis of Quadrature-Phase (or Imaginary) toward changes in MS.
% SENS_IP_con_cumu - Cumulative sensitivity analysis of In-Phase (or Real) toward changes in EC. 
% SENS_QP_con_cumu - Cumulative sensitivity analysis of Quadrature-Phase (or Imaginary) toward changes in EC.
% ErrorConIp - Error made on the sensitivity of In-Phase toward changes in EC.
% ErrorConQp - Error made on the sensitivity of Quadrature-Phase toward changes in EC.
% ErrorSusIp - Error made on the sensitivity of In-Phase toward changes in MS.
% ErrorSusQp - Error made on the sensitivity of Quadrature-Phase toward changes in MS.
%
%
% Functions:
%
% FDEM1DFWD_PM - Calculates the forward response (ppm) of a given layered 
%               half-space and loop-loop configuration. Calculation of the 
%               Hankel transform makes use of a digital filtering (Guptasarma 
%               and Singh filter). Use the Propagation Matrix approach.
% FDEM1DFWD_RC - Calculates the forward response (ppm) of a given layered 
%               half-space and loop-loop configuration. Calculation of the 
%               Hankel transform makes use of a digital filtering (Guptasarma 
%               and Singh filter). Use the Reflection Coefficient approach.
% FDEM1DSENS_PM - Calculates the sensitivity distribution of a given layered 
%                 soil medium and loop-loop configuration towards a certain 
%                 physical property using the brute-force or perturbation 
%                 method.
% FDEM1DSENS_RC - Calculates the sensitivity distribution of a given layered 
%                 soil medium and loop-loop configuration towards a certain 
%                 physical property using the brute-force or perturbation 
%                 method.
%
%
% [modifications]
%  OCT 2019 (JN) - For loop is replaced by arrayfun. manipulation of output
%                  matrix also in the Reflection Coefficient approach.
%  AUG 2019 (JN & LA) - For loop is replaced by arrayfun. manipulation of 
%                       output matrix so it work in the correlation.
%
% João Narciso and Leonardo Azevedo
% CERENA/IST, Portugal
% July 2019

%% ALLOCATE AND INTIALIZE
size_grid       = size(MSmodel);

%% FORWARD MODEL
if (strcmp(EMapprox ,'PM') == 1) 

    temp_i_i    = repmat((1:size_grid(1)), size(1:size_grid(2),2), 1);
    temp_j_j    = repmat((1:size_grid(2))', size(1:size_grid(1), 2),1);
    i           = temp_i_i(:);
    j           = temp_j_j(:);

    [SENS_IP_Sus, SENS_QP_Sus, SENS_IP_Con, SENS_QP_Con, FWD_IP, FWD_QP,...
        ErrorConIp, ErrorConQp, ErrorSusIp, ErrorSusQp] = arrayfun(@(i_, j_) FDEM1DSENS_PM_2SENS_array(i_,j_,ECmodel,MSmodel,S,M),i,j,'uniformoutput',0); % This has the wrong ordering

else

    temp_i_i    = repmat((1:size_grid(1)), size(1:size_grid(2),2), 1);
    temp_j_j    = repmat((1:size_grid(2))', size(1:size_grid(1), 2),1);
    i           = temp_i_i(:);
    j           = temp_j_j(:);

    [SENS_IP_Sus, SENS_QP_Sus, SENS_IP_Con, SENS_QP_Con, FWD_IP, FWD_QP,...
        ErrorConIp, ErrorConQp, ErrorSusIp, ErrorSusQp] = arrayfun(@(i_, j_) FDEM1DSENS_RC_2SENS_array(i_,j_,ECmodel,MSmodel,S,M),i,j,'uniformoutput',0); % This has the wrong ordering

%     for i = 1:size_grid(1)
%         for j = 1:size_grid(2)
% 
%            M.con = squeeze(ECmodel(i ,j , :))';  
%            M.sus = squeeze(MSmodel(i ,j , :))';              
% 
%            [FWD_IP(i,j,:), FWD_QP(i,j,:)] = FDEM1DFWD_RC(S,M);             % Reflection coefficient approach
% 
%            [SENS_IP_Sus(i, j, :), SENS_QP_Sus(i,j,:), ~] = FDEM1DSENS_RC(S,M,par); % Reflection coefficient approach 
%         end
%     end
end

FWD_IP  = cell2mat(FWD_IP)';
FWD_QP  = cell2mat(FWD_QP)';
ErrorConIp  = cell2mat(ErrorConIp)';
ErrorConQp  = cell2mat(ErrorConQp)';
ErrorSusIp  = cell2mat(ErrorSusIp)';
ErrorSusQp  = cell2mat(ErrorSusQp)';
SENS_IP_Sus = reshape(cell2mat(SENS_IP_Sus),size_grid(1), size_grid(2), size_grid(3));
% SENS_IP_Sus_cumu = reshape(cell2mat(SENS_IP_Sus_cumu),size_grid(1), size_grid(2), size_grid(3));
SENS_QP_Sus = reshape(cell2mat(SENS_QP_Sus),size_grid(1), size_grid(2), size_grid(3));
% SENS_QP_Sus_cumu = reshape(cell2mat(SENS_QP_Sus_cumu),size_grid(1), size_grid(2), size_grid(3));
SENS_IP_Con = reshape(cell2mat(SENS_IP_Con),size_grid(1), size_grid(2), size_grid(3));
% SENS_IP_Con_cumu = reshape(cell2mat(SENS_IP_Con_cumu),size_grid(1), size_grid(2), size_grid(3));
SENS_QP_Con = reshape(cell2mat(SENS_QP_Con),size_grid(1), size_grid(2), size_grid(3));
% SENS_QP_Con_cumu = reshape(cell2mat(SENS_QP_Con_cumu),size_grid(1), size_grid(2), size_grid(3));


%% Make figures of Sensitivity Analysis

% Figure of Sensitivity Analysis for IP and Electrical Conductivity
% 
% figure(1);
% Threshold = 0.7;
% subplot(1,2,1);
% Sens = (squeeze(SENS_IP_Con(1,:,1:end-1)))';
% SensCumu = (squeeze(SENS_IP_Con_cumu(1,:,1:end-1)))';
% idx = SensCumu <= 0.70;
% Sens70_Use = Sens.*idx;
% Sens70_NoUse = Sens.*~idx;
% Sens70_Use(Sens70_Use==0) = nan;
% Sens70_NoUse(Sens70_NoUse==0) = nan;
% SensCumu70Top = SensCumu.*idx;
% SensCumu70Base = SensCumu.*~idx;
% SensCumu70Top(SensCumu70Top==0) = nan;
% SensCumu70Base(SensCumu70Base==0) = nan;
% hold on;
%     plot(Sens70_Use,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(Sens70_NoUse,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     legend(['Mean Error < ',num2str(abs(mean(ErrorConIp)))], 'Location', 'southeast'); legend boxoff;
%     xlabel('Sensitivity'); ylabel('Depth (m)'); 
% hold off;
% subplot(1,2,2);
% hold on;
%     plot(SensCumu70Top,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(SensCumu70Base,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     plot([Threshold Threshold],[0 6],'r-', 'LineWidth', 2); 
%     text(0.7,5,'\rightarrow Threshold');
%     text(0.15,5.5,['Mean Error < ',num2str(abs(mean(ErrorConIp)))]);
%     xlabel('Cumulative Sensitivity');  
% hold off;
% text(-1, -0.3,['Sensitivity Analysis of IP for Electrical Conductivity.  Offset: ',num2str(S.x),' m'], 'FontSize', 16);
% savefig(figure(1),[pwd '/output/images/SensIP',num2str(par.ParSens1),'_',num2str(S.ori),'_Off',num2str(S.x),'_It',num2str(It),'_Sim',num2str(Sim),'.fig'],'compact');
% close all
% 
% 
% % Figure of Sensitivity Analysis for QP and Electrical Conductivity
% 
% figure(2); 
% Threshold = 0.7;
% subplot(1,2,1);
% Sens = (squeeze(SENS_QP_Con(1,:,1:end-1)))';
% SensCumu = (squeeze(SENS_QP_Con_cumu(1,:,1:end-1)))';
% idx = SensCumu <= 0.70;
% Sens70_Use = Sens.*idx;
% Sens70_NoUse = Sens.*~idx;
% Sens70_Use(Sens70_Use==0) = nan;
% Sens70_NoUse(Sens70_NoUse==0) = nan;
% SensCumu70Top = SensCumu.*idx;
% SensCumu70Base = SensCumu.*~idx;
% SensCumu70Top(SensCumu70Top==0) = nan;
% SensCumu70Base(SensCumu70Base==0) = nan;
% hold on;
%     plot(Sens70_Use,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(Sens70_NoUse,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     legend(['Mean Error < ',num2str(abs(mean(ErrorConQp)))], 'Location', 'southeast'); legend boxoff;
%     xlabel('Sensitivity'); ylabel('Depth (m)'); 
% hold off;
% subplot(1,2,2);
% hold on;
%     plot(SensCumu70Top,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(SensCumu70Base,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     plot([Threshold Threshold],[0 6],'r-', 'LineWidth', 2); 
%     text(0.7,5,'\rightarrow Threshold');
%     text(0.15,5.5,['Mean Error < ',num2str(abs(mean(ErrorConQp)))]);
%     xlabel('Cumulative Sensitivity');  
% hold off;
% text(-1, -0.3,['Sensitivity Analysis of QP for Electrical Conductivity.  Offset: ',num2str(S.x),' m'], 'FontSize', 16);
% savefig(figure(2),[pwd '/output/images/SensQP',num2str(par.ParSens1),'_',num2str(S.ori),'_Off',num2str(S.x),'_It',num2str(It),'_Sim',num2str(Sim),'.fig'],'compact');
% close all
% 
% % Figure of Sensitivity Analysis for IP and Magnetic Susceptibility
% 
% figure(3); 
% Threshold = 0.7;
% subplot(1,2,1);
% Sens = (squeeze(SENS_IP_Sus(1,:,1:end-1)))';
% SensCumu = (squeeze(SENS_IP_Sus_cumu(1,:,1:end-1)))';
% idx = SensCumu <= 0.70;
% Sens70_Use = Sens.*idx;
% Sens70_NoUse = Sens.*~idx;
% Sens70_Use(Sens70_Use==0) = nan;
% Sens70_NoUse(Sens70_NoUse==0) = nan;
% SensCumu70Top = SensCumu.*idx;
% SensCumu70Base = SensCumu.*~idx;
% SensCumu70Top(SensCumu70Top==0) = nan;
% SensCumu70Base(SensCumu70Base==0) = nan;
% hold on;
%     plot(Sens70_Use,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(Sens70_NoUse,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     legend(['Mean Error < ',num2str(abs(mean(ErrorSusIp)))], 'Location', 'southeast'); legend boxoff;
%     xlabel('Sensitivity'); ylabel('Depth (m)'); 
% hold off;
% subplot(1,2,2);
% hold on;
%     plot(SensCumu70Top,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(SensCumu70Base,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     plot([Threshold Threshold],[0 6],'r-', 'LineWidth', 2); 
%     text(0.7,5,'\rightarrow Threshold');
%     text(0.15,5.5,['Mean Error < ',num2str(abs(mean(ErrorSusIp)))]);
%     xlabel('Cumulative Sensitivity');  
% hold off;
% text(-1, -0.3,['Sensitivity Analysis of IP for Magnetic Susceptibility.  Offset: ',num2str(S.x),' m'], 'FontSize', 16);
% savefig(figure(3),[pwd '/output/images/SensIP',num2str(par.ParSens2),'_',num2str(S.ori),'_Off',num2str(S.x),'_It',num2str(It),'_Sim',num2str(Sim),'.fig'],'compact');
% close all
% 
% 
% % Figure of Sensitivity Analysis for QP and Magnetic Susceptibility
% 
% figure(4);
% Threshold = 0.7;
% Sens = (squeeze(SENS_QP_Sus(1,:,1:end-1)))';
% SensCumu = (squeeze(SENS_QP_Sus_cumu(1,:,1:end-1)))';
% idx = SensCumu <= 0.70;
% Sens70_Use = Sens.*idx;
% Sens70_NoUse = Sens.*~idx;
% Sens70_Use(Sens70_Use==0) = nan;
% Sens70_NoUse(Sens70_NoUse==0) = nan;
% SensCumu70Top = SensCumu.*idx;
% SensCumu70Base = SensCumu.*~idx;
% SensCumu70Top(SensCumu70Top==0) = nan;
% SensCumu70Base(SensCumu70Base==0) = nan;
% subplot(1,2,1);
% hold on;
%     plot(Sens70_Use,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(Sens70_NoUse,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     legend(['Mean Error < ',num2str(abs(mean(ErrorSusQp)))], 'Location', 'southeast'); legend boxoff;
%     xlabel('Sensitivity'); ylabel('Depth (m)'); 
% hold off;
% subplot(1,2,2);
% hold on;
%     plot(SensCumu70Top,cumsum(M.thick(1:end-1)),'g.'); set(gca, 'YDir','reverse');
%     plot(SensCumu70Base,cumsum(M.thick(1:end-1)),'k.'); set(gca, 'YDir','reverse');
%     plot([Threshold Threshold],[0 6],'r-', 'LineWidth', 2); 
%     text(0.7,5,'\rightarrow Threshold');
%     text(0.15,5.5,['Mean Error < ',num2str(abs(mean(ErrorSusQp)))]);
%     xlabel('Cumulative Sensitivity');  
% hold off;
% text(-1, -0.3,['Sensitivity Analysis of QP for Magnetic Susceptibility.  Offset: ',num2str(S.x),' m'], 'FontSize', 16);
% savefig(figure(4),[pwd '/output/images/SensQP',num2str(par.ParSens2),'_',num2str(S.ori),'_Off',num2str(S.x),'_It',num2str(It),'_Sim',num2str(Sim),'.fig'],'compact');
% close all
% 
end


