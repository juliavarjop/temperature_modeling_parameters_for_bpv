function [TmodM,TmodS,RMSdifference] = calculateTmodsAtConditions_VBPV_v1(convectionCoeffs,conditions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tsimulated calculation %%%

errorCountLimit=5;

%%% Parameters
% Natural constants
consts = ([]);
consts.('h') = 6.626*10^-34;    % Js
consts.('h_eV') = 4.136*10^-15; % eVs
consts.('c') = 2.998*10^8;      % m/s
consts.('kB') = 1.381*10^-23;   % J/K
consts.('q') = 1.602*10^-19;    % C
consts.('SB') = 5.670*10^-8;    % W/(m^2K^4);
% 
% Band gap
% Si
EgRef = 1.206; % eV, https://doi.org/10.1063/1.345414
EgRefT = 0; % K, https://doi.org/10.1063/1.345414
EgShiftVsT = -0.273*1e-3; % eV/K, https://doi.org/10.1063/1.345414

% Total efficiencies
EffRefT = 298.15; % K
EffRefG = 1000; % W/m^2
VocRefG = 0.6; % V
n=1.3;

% %%%% CBPV parameters 
% effRefFront = 0.195;
% effRefBack = 0.195*0.7;
% EffShiftVsTSi = -0.354; % %/K

% %%%% VBPV parameters 
effRefFront = 0.177;
effRefBack = 0.159;
EffShiftVsTSi = -0.4139; % %/K

% Spectral absorption coeff above and below the respective band gaps
absCoeffAboveGap = 0.95;
absCoeffBelowGap = 0.2;

% Irradiation data
AM15 = load('AM15.mat').AM15;

% Modelled wavelengths
dLambda = 1; % nm
lambda = 300:dLambda:4000; 
lambda = transpose(lambda);

% Modelled spectrum
lambda2 = AM15{:,1};                        % Wavelength (nm)
B = AM15{:,2};                              % Radiance (Wm^-2nm^-1)
Bsolar = spline(lambda2,B,lambda);
BsolarTot = sum(Bsolar*dLambda);

% Model validation with TUAS VBPV data %

% Initial temperature 
Tstc = 293.15; % K

% Folder
pathFolder ='';


% % Error log
% FolderErrorLog = [pathFolder,'\errorLog'];
% mkdir(FolderErrorLog);

% Weather conditions

I_Front_R=conditions.POA_front;
I_Rear_R=conditions.POA_rear;
v_wind=conditions.v_wind;
Tamb=conditions.Tamb;

% Temperature convergence
TendCriterion = 0.1; % K
countTlimit = 10;

a_coeff=convectionCoeffs(1);
b_coeff=convectionCoeffs(2);

disp(['a = ', num2str(convectionCoeffs(1)), ', b = ', num2str(convectionCoeffs(2))])


%%%% h = a_coeff * v + 2.8, where k_coeff = 1 %%%%

[m, ~] = size(conditions(:,:));

TaveCells_vs_Si = zeros(length(m),countTlimit);
TcellAveTop_vs_Si = zeros(length(m),countTlimit);
TcellAveBack_vs_Si = zeros(length(m),countTlimit);
TaveMod_vs_Si = zeros(length(m),countTlimit);
TmaxCells_vs_Si = zeros(length(m),countTlimit);

CutLineCellsResults_vs_Si = cell(length(m),countTlimit);
CutLineTopResults_vs_Si = cell(length(m),countTlimit);
CutLineBackResults_vs_Si = cell(length(m),countTlimit);
CutLineCellsXResults_vs_Si = cell(length(m),countTlimit);
CutLineTopXResults_vs_Si = cell(length(m),countTlimit);
CutLineBackXResults_vs_Si = cell(length(m),countTlimit);

% Loop over environmental conditions
for idxConditions = 1:m
    % Start time
    tStartRound = tic;

    %%% Loop over temperature iteration
    countT = 0;
    Tchange = 999;
    errorCount = 0;
    Tcell = Tstc; % K

    while (Tchange > TendCriterion && countT < countTlimit)
        try

            % Band gap at a temperature
            EgAtT = EgRef+EgShiftVsT*(Tcell-EgRefT);
            EgWl = consts.h_eV*consts.c*10^9/EgAtT; % nm

            % Efficiency at current temperature

            % Absorption
            A = ones(size(Bsolar))*absCoeffBelowGap;
            A(lambda<EgWl) = absCoeffAboveGap;

            % Total absorption coefficient
            % Spectrum weighted average of normalized wavelength specific absorption
            totAbsCoeff = sum(Bsolar.*A*dLambda)/BsolarTot;

            %%% Qcell Front
            % Electrical efficiency
            EffRefAtT = effRefFront*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((I_Front_R(idxConditions))/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);
            effEref = EffRefAtT/totAbsCoeff;

            % Heat production
            QprodCoeff = (1-effEref)*totAbsCoeff;

            % Irradiance effect
            Qcell_Front = I_Front_R(idxConditions)*QprodCoeff; % W/m^2

            %%% Qcell Back
            % Electrical efficiency
            EffRefAtT = effRefBack*(1+(n*consts.kB*EffRefT)/(VocRefG*consts.q)*log((I_Rear_R(idxConditions))/(EffRefG)))*(1+EffShiftVsTSi*(Tcell-EffRefT)/100);
            effEref = EffRefAtT/totAbsCoeff;

            % Heat production
            QprodCoeff = (1-effEref)*totAbsCoeff;

            % Irradiance effect
            Qcell_Back = I_Rear_R(idxConditions)*QprodCoeff; % W/m^2

            %%% Qcell Total
            Qcell=Qcell_Front+Qcell_Back;

            % Build Comsol model
            % VBPV panel COMSOL model
            [~,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,~,~,~,~,~,~,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = VBPV_SipanelTcomsolModel(Tamb(idxConditions),v_wind(idxConditions),Qcell,pathFolder,Tcell,a_coeff,b_coeff);
            
            % CBPV panel COMSOL model
            % [~,TcellAve,TcellAveTop,TcellAveBack,TmodAve,TcellMax,~,~,~,~,~,~,CutLineCells,CutLineTop,CutLineBack,CutLineCellsX,CutLineTopX,CutLineBackX] = MPV_CBPV_SipanelTcomsolModel(Tamb(idxConditions),v_wind(idxConditions),Qcell,pathFolder,Tcell,a_coeff,b_coeff);
            % Extract results
            TcellUpdate = TcellAve;

            % Temperature convergence?
            Tchange = abs(TcellUpdate-Tcell);

            % Update cell temperature
            Tcell = TcellUpdate;  % Cell average

            % T iteration limit reached
            countT = countT+1;
            if countT >= countTlimit
                warning('Tcell iteration limit reached.') % TO DO: save into log file
            end

            % Save results
            TaveCells_vs_Si(idxConditions,countT) = Tcell;
            TcellAveTop_vs_Si(idxConditions,countT) = TcellAveTop;
            TcellAveBack_vs_Si(idxConditions,countT) = TcellAveBack;
            TaveMod_vs_Si(idxConditions,countT) = TmodAve;
            TmaxCells_vs_Si(idxConditions,countT) = TcellMax;

            CutLineCellsResults_vs_Si{idxConditions,countT} = CutLineCells;
            CutLineTopResults_vs_Si{idxConditions,countT} = CutLineTop;
            CutLineBackResults_vs_Si{idxConditions,countT} = CutLineBack;
            CutLineCellsXResults_vs_Si{idxConditions,countT} = CutLineCellsX;
            CutLineTopXResults_vs_Si{idxConditions,countT} = CutLineTopX;
            CutLineBackXResults_vs_Si{idxConditions,countT} = CutLineBackX;


            errorCount = 0;
        catch ME
            errorCount = errorCount+1;

            % Write error to a file
            % fid = fopen([FolderErrorLog,'\errorCount_',num2str(errorCount),'.txt'],'w');
            % fprintf(fid,'%s',ME.getReport('extended','hyperlinks','off'));
            % fclose(fid);

            pause(5)
            if errorCount>errorCountLimit
                keyboard
            end
        end
    end
    % Progress report
    tElapsedRound = toc(tStartRound);
    disp(['Runtime, conditions (',num2str(idxConditions),'/',num2str(m),...
        '): ',num2str(tElapsedRound/60),' min.'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tsimulated vs Tmeasured error metrics %%%

%%% Find the last iteration and save the operating temperature as celsius degreees %%%

indDataPoint = 1;
for indm = 1:m 
    ind = find(TaveCells_vs_Si(indm,:),1,'last');

    indT = find(~cellfun('isempty', CutLineBackResults_vs_Si(indm,:)), 1, 'last'); 
    indx = find(~cellfun('isempty', CutLineBackXResults_vs_Si(indm,:)), 1, 'last'); 

    x=CutLineBackXResults_vs_Si{indm,indx}*100;
    y=CutLineBackResults_vs_Si{indm,indT};

    r = 1;
    k = 1;
    limitMIN = 11;
    limitMAX = 12;
    listMIN = [];
    listMAX = [];

    for i = 1:limitMIN
        idx = find(x>=r*(1+15.6) & x<=(r+1)*1+r*15.6);
        yidx = y(idx);
        r = r + 1;

        listMIN = [listMIN;{yidx}];
    end

    for j = 1:limitMAX
        idx2 = find(x>=k*1+(k-1)*15.6 & x<=k*(1+15.6));
        yidx = y(idx2);
        k = k + 1;

        listMAX = [listMAX;{yidx}];
    end

    % Store results 
    
    % The maximum temperature of the module back surface in varying combinations of ambient conditions
    TmaxBack(indDataPoint,1) = I_Front_R(indm);
    TmaxBack(indDataPoint,2) = I_Rear_R(indm);
    TmaxBack(indDataPoint,3) = Tamb(indm);
    TmaxBack(indDataPoint,4) = v_wind(indm);
    TmaxBack(indDataPoint,5)=max([listMAX{:}])-273.15;

    % The minimum temperature of the module back surface in varying combinations of ambient conditions
    TminBack(indDataPoint,1) = I_Front_R(indm);
    TminBack(indDataPoint,2) = I_Rear_R(indm);
    TminBack(indDataPoint,3) = Tamb(indm);
    TminBack(indDataPoint,4) = v_wind(indm);
    TminBack(indDataPoint,5) = min([listMIN{:}])-273.15;
    % 
    % The average temperature of the cells in varying combinations of ambient conditions
    TaveCells(indDataPoint,1) = I_Front_R(indm);
    TaveCells(indDataPoint,2) = I_Rear_R(indm);
    TaveCells(indDataPoint,3) = Tamb(indm);
    TaveCells(indDataPoint,4) = v_wind(indm);
    TaveCells(indDataPoint,5) = TaveCells_vs_Si(indm,ind)-273.15;

    % The average temperature of the back surface in varying combinations of ambient conditions
    TaveBack(indDataPoint,1) = I_Front_R(indm);
    TaveBack(indDataPoint,2) = I_Rear_R(indm);
    TaveBack(indDataPoint,3) = Tamb(indm);
    TaveBack(indDataPoint,4) = v_wind(indm);
    TaveBack(indDataPoint,5) = TcellAveBack_vs_Si(indm,ind)-273.15;
    
    indDataPoint = indDataPoint+1;

end

%%% Calculate RMSE %%%

TmodM=conditions.T2(:,1);
TmodS=TminBack(:,5);

filt=isfinite(TmodM) & isfinite(TmodS);

diff = TmodM(filt) - TmodS(filt);
RMSdifference = sqrt(mean(diff.^2));

%disp(['RMSE = ', num2str(RMSdifference)])


end