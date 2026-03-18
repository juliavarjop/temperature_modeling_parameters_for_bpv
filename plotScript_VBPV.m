% Plot script to simulate the operating temperature of VBPV panel
% in given weather conditions using optimized (or unoptimized) convective heat transfer coefficients

% Clear workspace
% clear all;

% Initial data
dataInitial=load('VBPV_TUAS_2018to2020_measdata_and_simirradiancedata_v2.mat').VBPV_TUAS_2018to2020_measdata_and_simirradiancedata_v2;

Ta = 5:5:30; 
Irradiance = 200:100:1000; 
WS = 1:1:5; 

start_WD=0;
end_WD=360;

fitCoeffs = zeros(length(start_WD));

fitCoeffs(:,1)=1.1113; % fitted (optimized) coefficient a
fitCoeffs(:,2)=2.9940; % fitted (optimized) coefficient b

%a=1.7669, b=7.9677 h coeffs of MPV panel (with snow filter)
%a=2.4147, b=2.9087 h coeffs of CBPV panel (with snow filter)
%a=1.1113, b=2.9940 h coeffs of VBPV panel (with snow filter)

fitVariables0(1,1)=3; % initial (unoptimized) coefficient a
fitVariables0(1,2)=2.8; % initial (unoptimized) coefficient b

pathFolder ='';

TmodSfit_table = cell(length(start_WD)); % simulated module temperature, fitted (optimized) coefficients
TmodMfit_table = cell(length(start_WD));
RMSdifferenceFit_table = cell(length(start_WD));
TmodS0_table = cell(length(start_WD));  % simulated module temperature, initial (unoptimized) coefficients
TmodM0_table = cell(length(start_WD));
RMSdifference0_table = cell(length(start_WD));

for ind_WD=1:length(start_WD)

    % Filter data
    WD_filter = dataInitial.windDirection>=start_WD(ind_WD) & dataInitial.windDirection<=end_WD(ind_WD);
    data=dataInitial(WD_filter,:);

    TambExperimentalData = data.TEMP;
    POAirradianceData = data.I_Rear_R+data.I_Front_R;
    POAWestirradianceData = data.I_Rear_R;
    POAEastirradianceData = data.I_Front_R;
    wsExperimentalData = data.windSpeed;
    TexperimentalData = data.thermocouple2;
    WindDirectionExperimentalData = data.windDirection;

    % Extract data
    DataConditions = zeros(length(WS)*length(Ta)*length(Irradiance),3);
    ConditionsClosestTambExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestPOAtot = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestPOAWest = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestPOAEast = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestWsExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    validationConditionsClosestWDExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);

    DataTexpClosest = zeros(length(WS)*length(Ta)*length(Irradiance),1);

    indDataPoint = 1;
    for indWS = 1:length(WS)
        for indIrr = 1:length(Irradiance)
            for indTamb = 1:length(Ta)

                % Store conditions
                DataConditions(indDataPoint,1) = Ta(indTamb);
                DataConditions(indDataPoint,2) = Irradiance(indIrr);
                DataConditions(indDataPoint,3) = WS(indWS);

                % Create validation set
                conditionsDiff = ...
                    [TambExperimentalData,POAirradianceData,wsExperimentalData]-...
                    [Ta(indTamb),Irradiance(indIrr),WS(indWS)];
                [~,indClosest] = min(vecnorm(conditionsDiff,2,2));
                ConditionsClosestTambExp(indDataPoint) = TambExperimentalData(indClosest);
                ConditionsClosestPOAtot(indDataPoint) = POAirradianceData(indClosest);
                ConditionsClosestPOAWest(indDataPoint) = POAWestirradianceData(indClosest);
                ConditionsClosestPOAEast(indDataPoint) = POAEastirradianceData(indClosest);
                ConditionsClosestWsExp(indDataPoint) = wsExperimentalData(indClosest);
                validationConditionsClosestWDExp(indDataPoint) = WindDirectionExperimentalData(indClosest);
                DataTexpClosest(indDataPoint) = TexperimentalData(indClosest);
                indDataPoint = indDataPoint+1;
            end
        end
    end


    VBPV_InputData=table(ConditionsClosestPOAtot,...
        ConditionsClosestPOAWest,ConditionsClosestPOAEast,...
        ConditionsClosestWsExp,validationConditionsClosestWDExp,...
        ConditionsClosestTambExp,DataTexpClosest,...
        'VariableNames',["POA","POA_rear","POA_front","v_wind","v_direct","Tamb","T2"]);

    [~,ib] = unique(VBPV_InputData);
    VBPV_InputData = VBPV_InputData(ib,:);
    
    % Get data with the fitted (optimized) and intial (unoptimized) coefficients
      
    [TmodMfit,TmodSfit,RMSdifferenceFit] = calculateTmodsAtConditions_VBPV_v1(fitCoeffs(ind_WD,1:2),VBPV_InputData);
    [TmodM0,TmodS0,RMSdifference0] = calculateTmodsAtConditions_VBPV_v1(fitVariables0,VBPV_InputData);

    TmodSfit_table{ind_WD,1} = TmodSfit; % simulated module temperature, fitted (optimized) coefficients
    TmodMfit_table{ind_WD,1} = TmodMfit;
    RMSdifferenceFit_table{ind_WD,1} = RMSdifferenceFit;
    TmodS0_table{ind_WD,1} = TmodS0; % simulated module temperature, initial (unoptimized) coefficients
    TmodM0_table{ind_WD,1} = TmodM0;
    RMSdifference0_table{ind_WD,1} = RMSdifference0;

end

VBPV_TminBack_fitCoeffs=TmodSfit_table{1,1};

% Save workspace
save([pathFolder,'\workspace.mat'])
