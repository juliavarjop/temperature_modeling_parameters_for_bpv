
% Script to find the optimized convective heat transfer coefficient for the VBPV panel

% Folder
pathFolder ='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Data selection %%%
dataInitial=load('VBPV_TUAS_2018to2020_measdata_and_simirradiancedata_v2.mat').VBPV_TUAS_2018to2020_measdata_and_simirradiancedata_v2;

% Conditions
Ta = 5:5:30; 
Irradiance = 200:100:1000;
WS = 1:1:5; 

start_WD=0; %[0 45 90 135 180 225 270 315];
end_WD=360; %start_WD+45;

parameter_a = zeros(length(start_WD));
parameter_b = zeros(length(start_WD));

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
                ConditionsClosestPOAWest(indDataPoint) = POAWestirradianceData(indClosest);
                ConditionsClosestPOAEast(indDataPoint) = POAEastirradianceData(indClosest);
                ConditionsClosestWsExp(indDataPoint) = wsExperimentalData(indClosest);
                ConditionsClosestWDExp(indDataPoint) = WindDirectionExperimentalData(indClosest);
                DataTexpClosest(indDataPoint) = TexperimentalData(indClosest);
                indDataPoint = indDataPoint+1;
            end
        end
    end


    VBPV_InputData_10min_v1=table(ConditionsClosestPOAtot,...
        ConditionsClosestPOAWest,ConditionsClosestPOAEast,...
        ConditionsClosestWsExp,ConditionsClosestWDExp,...
        ConditionsClosestTambExp,DataTexpClosest,...
        'VariableNames',["POA","POAWest","POAEast","v_wind","v_direct","Tamb","T2"]);

    [~,ib] = unique(VBPV_InputData_10min_v1);
    InputData = VBPV_InputData_10min_v1(ib,:);

    % Fit
    % Initial guesses (example; in reality, make as good guesses as possible to
    % limit the number of iterations and therefore compuation time)
    a0 = 1;
    b0 = 3;
    fitVariables0 = [a0,b0];
    % Constrains (if any? ask about formatting if needed)
    AConstr = [];
    bConstr = [];
    Aeq = [];
    beq = [];
    % Bounds (for variables a,b): 0 <= functionVariables <= 10
    lb = [0,0];
    ub = [10,10];
    % Non linear conditions
    nonlcon = [];
    % Optimization options: set ending criteria
    % Note:
    % - Suitable FunctionTolerance depends on the chosen metric (here, RMS) and
    %   how close we need to get
    % - StepTolerance referes to the changes in a and b
    % - 'Display','iter' allows to follow how the optimization is progressing
    %   by each iteration
    options = optimoptions(@fmincon,'Display','iter-detailed','StepTolerance',2e-2,...
        'FunctionTolerance',5e-2,'DiffMinChange', 1e-2);
    % options = optimoptions(@fmincon,'Display','iter');

    % Fit
    % fitCoeffs = fmincon(@(fitVariables) ...
    %     calculateRMSdifference_v1(fitVariables,InputData),...
    %     fitVariables0,AConstr,bConstr,Aeq,beq,lb,ub,nonlcon,options);

    [fitCoeffs,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(fitVariables) ...
        calculateRMSdifference_VBPV_v1(fitVariables,InputData),...
        fitVariables0,AConstr,bConstr,Aeq,beq,lb,ub,nonlcon,options);


    % Print fitted values
    disp(['Fit a: ', num2str(fitCoeffs(1)),', fit b: ', num2str(fitCoeffs(2))])

    parameter_a(ind_WD,1) = fitCoeffs(1);
    parameter_b(ind_WD,1) = fitCoeffs(2);

    % h = parameter_a*v + parameter_b; 

end

% Save workspace
save([pathFolder,'\workspace.mat'])
