folder='';

%%% INPUT DATA %%%

data=load('MPV_InputData.mat').MPV_InputData;

Tamb=data.Tamb;
POA=data.POA;
WS=data.v_wind;

T_measured=load('MPV_TmaxBack_fitCoeffs.mat').MPV_TmaxBack_fitCoeffs;


%%% SANDIA MODEL PARAMETERS %%%

y=log((T_measured-Tamb)./POA);

x1 = WS>=0;
x2 = WS<=4;
x3 = WS>=4;
x4 = WS<=10;

% fit between 0-4 m/s
Sandia_p_0to4 = polyfit(WS(x1 & x2),y(x1 & x2),1);

% fit between 4-10 m/s 
Sandia_p_4to10 = polyfit(WS(x3 & x4),y(x3 & x4),1);

% fit between 0-10 m/s
Sandia_p_0to10 = polyfit(WS(x1 & x4),y(x1 & x4),1)

a=Sandia_p_0to10(:,2);
b=Sandia_p_0to10(:,1);

%%% FIGURE %%%

width = 700; 
height = 660; 

figure('Position', [100, 100, width, height]);

% Line shows between 0-4 m/s
WS_0to4 = WS(WS>=0 & WS<=4);
y_0to4 = Sandia_p_0to4(:,1) * WS_0to4+ Sandia_p_0to4(:,2);

% Line shows between 4-10 m/s
WS_4to10 = WS(WS>=4 & WS<=10);
y_4to10 = Sandia_p_4to10(:,1) * WS_4to10+ Sandia_p_4to10(:,2);

% Line shows between 0-10 m/s
WS_0to10 = WS(WS>=0 & WS<=10);
y_0to10 = Sandia_p_0to10(:,1) * WS_0to10+ Sandia_p_0to10(:,2);

plot8=plot(WS,y,'.','Color',[0, 0.5, 0.9],'MarkerSize',25,'LineWidth',1.25);
hold on
plot9=plot(WS_0to10,y_0to10,'Color',[0 0 0],'LineWidth',1.25);

% limits
ylim([-4.7 -3.1])
xlim([0 12])
yticks(-4.7:0.4:-3.1)
xticks(0:3:12)
set(gca, 'FontSize', 22); 

hold off
set(gca, 'LineWidth', 1);
grid on;

xlabel('Wind speed (m/s)', FontSize=22);
ylabel('log((T_m-T_a_m_b)/G)', FontSize=22);
hold on;

ax1 = gca;

lgd1=legend(ax1, [plot9], ['0-10 m/s: y=',num2str(round(Sandia_p_0to10(1,1),4)),'x',num2str(round(Sandia_p_0to10(1,2),3))], 'FontSize', 22, 'Location', 'southwest');

hold off;

%%% SAVE FIGURE %%% 
% saveas(gcf,[folder,'\sandia_simulated_params_fmi_v2'],'epsc')
% saveas(gcf,[folder,'\\sandia_simulated_params_fmi_v2.png'])
% saveas(gcf,[folder,'\\\sandia_simulated_params_fmi_v2.fig'])


%%% FAIMAN MODEL PARAMETERS %%%

y = POA./(T_measured-Tamb);

x1 = WS>=0;
x4 = WS<=10;

% fit between 0-10 m/s
Faiman_p_0to10 = polyfit(WS(x1 & x4),y(x1 & x4),1)

U_L1=Faiman_p_0to10(:,1);
U_L0=Faiman_p_0to10(:,2);

%%% FIGURE %%%

width = 700; 
height = 660; 

figure('Position', [100, 100, width, height]);

% Line shows between 0-10 m/s
WS_0to10 = WS(WS>=0 & WS<=10);
y_0to10 = Faiman_p_0to10(:,1) * WS_0to10+ Faiman_p_0to10(:,2);

% lines
plot8=plot(WS,y,'.','Color',[0, 0.5, 0.9],'MarkerSize',25,'LineWidth',1.25);
hold on
plot9=plot(WS_0to10,y_0to10,'Color',[0 0 0],'LineWidth',1.25);

% % limits
xlim([0 12])
ylim([20 80])
xticks(0:3:12)
yticks(20:15:100)
set(gca, 'FontSize', 22); 

hold off
set(gca, 'LineWidth', 1);
grid on;

xlabel('Wind speed (m/s)', FontSize=22);
ylabel('G/(T_m-T_a_m_b)', FontSize=22);
hold on;

ax1 = gca;

lgd1=legend(ax1, [plot9], ['0-10 m/s: y=',num2str(round(Faiman_p_0to10(1,1),4)),'x+',num2str(round(Faiman_p_0to10(1,2),3))], 'FontSize', 22, 'Location', 'southeast');

hold off;

%%% SAVE FIGURE %%% 
% saveas(gcf,[folder,'\faiman_simulated_params_fmi_v2'],'epsc')
% saveas(gcf,[folder,'\\faiman_simulated_params_fmi_v2.png'])
% saveas(gcf,[folder,'\\\faiman_simulated_params_fmi_v2.fig'])


%%%  PVSYST MODEL PARAMETERS %%%

at=0.9;
y = POA.*at*(1-0.1491)./(T_measured-Tamb); %

x1 = WS>=0;
x4 = WS<=10;

% fit between 0-10 m/s
PVsyst_p_0to10 = polyfit(WS(x1 & x4),y(x1 & x4),1)

U1=PVsyst_p_0to10(:,1);
U0=PVsyst_p_0to10(:,2);

%%% FIGURE %%%

width = 700; 
height = 660; 

figure('Position', [100, 100, width, height]);

% Line shows between 0-10 m/s
WS_0to10 = WS(WS>=0 & WS<=10);
y_0to10 = PVsyst_p_0to10(:,1) * WS_0to10+ PVsyst_p_0to10(:,2);

% lines
plot8=plot(WS,y,'.','Color',[0, 0.5, 0.9],'MarkerSize',25,'LineWidth',1.25);
hold on
plot9=plot(WS_0to10,y_0to10,'Color',[0 0 0],'LineWidth',1.25);

% limits
xlim([0 12])
ylim([10 70])
xticks(0:3:12)
yticks(10:15:100)
set(gca, 'FontSize', 22); 

hold off
set(gca, 'LineWidth', 1);
grid on;

xlabel('Wind speed (m/s)', FontSize=22);
ylabel('^{[G \alpha_{tot} (1-PCE_{ref})]}/_{(T_{m}-T_{amb})}', FontSize=22);

hold on;

ax1 = gca;

lgd1=legend(ax1, [plot9], ['0-10 m/s: y=',num2str(round(PVsyst_p_0to10(1,1),4)),'x+',num2str(round(PVsyst_p_0to10(1,2),3))], 'FontSize', 22, 'Location', 'southwest');

hold off;

%%% SAVE FIGURE %%% 
% saveas(gcf,[folder,'\pvsyst_simulated_params_fmi_v2'],'epsc')
% saveas(gcf,[folder,'\\pvsyst_simulated_params_fmi_v2.png'])
% saveas(gcf,[folder,'\\\pvsyst_simulated_params_fmi_v2.fig'])