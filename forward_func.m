function d_pred =  forward_func(ECModel, MSModel)
% EM Approx. RC (Reflection Coefficients) or PM (Propagation Matrix)
EMApprox = 'RC';       % RC ou PM

% SENSOR CHARACTERISTICS (S structure)
% Model characteristics (M structure)

LayerThickness = 0.10;
nx = 400; nz = 40;

SStruct.x1 = [1 2]; % x-coordinate receiver (m) for HCP
SStruct.x2 = [1.1 2.1];  % x-coordinate receiver (m) for PRP
SStruct.y = 0; % y-coordinate receiver (m)
SStruct.z = -0.15; % z-coordinate receiver (m) - positive z-axis pointed down
SStruct.height = 0.15; % Height of transmitter (m)
SStruct.freq = 9000; % Frequency (Hz)
SStruct.mom = 1; % Transmitter moment (A.m^2)
SStruct.ori1 = 'ZZ'; % Coil orientation for x first line (Transmitter(X,Y,Z),Receiver(X,Y,Z))
SStruct.ori2 = 'ZX'; % Coil orientation for x second line (Transmitter(X,Y,Z),Receiver(X,Y,Z))

% Allocate Susceptibility, Conductivity, Permittivity and Layer thickness
MStruct.perm = ones(1, nz) .* 1e-12; % Permittivity of layer(s) (F/m)
MStruct.thick = repmat(LayerThickness, 1, nz); % Layer(s) thickness (m)

MSModel = reshape(MSModel, [1, nx, nz]);
ECModel = reshape(ECModel, [1, nx, nz]);
MStruct.sus = MSModel; % Susceptibility of layer(s) (-)
MStruct.con = ECModel;  % Conductivity of layer(s) (S/m)

%
offsetReal1 = [1 2];
offsetReal2 = [1.1 2.1];
OffsetTotal = [offsetReal1 offsetReal2];
noff = length(OffsetTotal);
ParSens = 'cons';
FWD_IP_off1 = zeros(length(offsetReal1), nx);
FWD_QP_off1 = zeros(length(offsetReal1), nx);

FWD_IP_off2 = zeros(length(offsetReal2), nx);
FWD_QP_off2 = zeros(length(offsetReal2), nx);

MSModel = MStruct.sus;
ECModel = MStruct.con;
if (strcmp(SStruct.ori1 ,'ZZ') == 1)  || (strcmp(SStruct.ori1 ,'YY') == 1)  || (strcmp(SStruct.ori1 ,'XX') == 1)
    SStruct.ori = SStruct.ori1;

    for noff1 = 1:2      
        SStruct.x = offsetReal1(1, noff1);
        [FWD_IP_off1(noff1, :), FWD_QP_off1(noff1, :)] = FWD_EM_2sens(SStruct, MStruct, ECModel, MSModel, ParSens, EMApprox);
    end
end


if (strcmp(SStruct.ori2 ,'ZX') == 1) || (strcmp(SStruct.ori2 ,'XZ') == 1) || (strcmp(SStruct.ori2 ,'XY') == 1)

    SStruct.ori = SStruct.ori2;

    for noff2 = 1:2

        SStruct.x = offsetReal2(1, noff2);

        [FWD_IP_off2(noff2, :), FWD_QP_off2(noff2, :)] = FWD_EM_2sens(SStruct, MStruct, ECModel, MSModel, ParSens, EMApprox);
    end
end

d_pred = [FWD_IP_off1; FWD_IP_off2; FWD_QP_off1; FWD_QP_off2];
d_pred = d_pred';
d_pred = d_pred(:);
