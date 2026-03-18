% Plot script to simulate the operating temperature of MPV panel
% in given weather conditions using optimized (or unoptimized) convective heat transfer coefficients

% Clear workspace
% clear all;

% Initial data
dataInitial=load('MPV_FMI_2015to2021_measdata_and_simirradiancedata.mat').MPV_FMI_2015to2021_measdata_and_simirradiancedata;

Ta = 20; %5:5:30; 
Irradiance = 400; %200:100:1000;
WS = 2; %1:1:10; 

start_WD=0;
end_WD=360;

fitCoeffs = zeros(length(start_WD));

fitCoeffs(:,1)=1.7669; % fitted (optimized) coefficient a
fitCoeffs(:,2)=7.9677; % fitted (optimized) coefficient b

%a=1.7669, b=7.9677 h coeffs of MPV panel (with snow filter)
%a=2.4147, b=2.9087 h coeffs of CBPV panel (with snow filter)
%a=1.1113, b=2.9940 h coeffs of VBPV panel (with snow filter)

fitVariables0(1,1)=3; % initial (unoptimized) coefficient a
fitVariables0(1,2)=2.8; % initial (unoptimized) coefficient b

pathFolder ='';

TmodSfit_table = cell(length(start_WD)); % simulated module temperature, fitted (optimized) coefficients
TmodMfit_table = cell(length(start_WD));
RMSdifferenceFit_table = cell(length(start_WD));
TmodS0_table = cell(length(start_WD)); % simulated module temperature, initial (unoptimized) coefficients
TmodM0_table = cell(length(start_WD));
RMSdifference0_table = cell(length(start_WD));


for ind_WD=1:length(start_WD)

    % Filter data
    WD_filter = dataInitial.WD_PT10M_AVG>=start_WD(ind_WD) & dataInitial.WD_PT10M_AVG<end_WD(ind_WD);
    data=dataInitial(WD_filter,:);

    TambExperimentalData = data.TA_PT1M_AVG31;
    POAirradianceData = data.I_Front_R;
    wsExperimentalData = data.WS_PT10M_AVG;
    TexperimentalData = data.TTECH_PT1M_AVG33; % SW
    WindDirectionExperimentalData = data.WD_PT10M_AVG;

    % Extract data
    DataConditions = zeros(length(WS)*length(Ta)*length(Irradiance),3);
    ConditionsClosestTambExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestPOAtot = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestWsExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);
    ConditionsClosestWDExp = zeros(length(WS)*length(Ta)*length(Irradiance),1);

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
                ConditionsClosestWsExp(indDataPoint) = wsExperimentalData(indClosest);
                ConditionsClosestWDExp(indDataPoint) = WindDirectionExperimentalData(indClosest);
                DataTexpClosest(indDataPoint) = TexperimentalData(indClosest);
                indDataPoint = indDataPoint+1;
            end
        end
    end


    MPV_InputData=table(ConditionsClosestPOAtot,...
        ConditionsClosestWsExp,ConditionsClosestWDExp,...
        ConditionsClosestTambExp,DataTexpClosest,...
        'VariableNames',["POA","v_wind","v_direct","Tamb","T"]);

    [~,ib] = unique(MPV_InputData);
    MPV_InputData = MPV_InputData(ib,:);

    
    % Get data with the fitted (optimized) and intial (unoptimized) coefficients
    
    [TmodMfit,TmodSfit,RMSdifferenceFit] = calculateTmodsAtConditions_MPV_v1(fitCoeffs(ind_WD,1:2),MPV_InputData);
    [TmodM0,TmodS0,RMSdifference0] = calculateTmodsAtConditions_MPV_v1(fitVariables0,MPV_InputData);

    TmodSfit_table{ind_WD,1} = TmodSfit; % simulated module temperature, fitted (optimized) coefficients
    TmodMfit_table{ind_WD,1} = TmodMfit;
    RMSdifferenceFit_table{ind_WD,1} = RMSdifferenceFit;
    TmodS0_table{ind_WD,1} = TmodS0; % simulated module temperature, initial (unoptimized) coefficients
    TmodM0_table{ind_WD,1} = TmodM0;
    RMSdifference0_table{ind_WD,1} = RMSdifference0;

end

MPV_TmaxBack_fitCoeffs=TmodSfit_table{1,1};

% Save workspace
save([pathFolder,'\workspace.mat'])
