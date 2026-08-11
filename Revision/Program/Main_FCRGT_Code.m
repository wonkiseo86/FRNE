% Main Code for Functional Cointegrating Regression between GRP Growth and Temperature Anomalies
%
% Final code for the JASA revision submission
% August 2026
% Kyungsik Nam (ksnam@hufs.ac.kr)

%% Set-up and Data Loading
clear; clc; close all;

code_dir = ['C:\Users\noran\Dropbox\Academic_Life(201808-Present)\' ...
    '2_Main_Research\0_Journal_Pub\2026_Submission\202602_Sector_CC\' ...
    '202509_P1_Code\Program\202608_JASA_Revision_Code'];
if ~isfolder(code_dir)
    code_file = mfilename('fullpath');
    if isempty(code_file) && usejava('desktop')
        code_file = matlab.desktop.editor.getActiveFilename;
    end
    if isempty(code_file)
        error(['The code folder could not be identified. Update code_dir in ' ...
            'the Set-up and Data Loading section.']);
    end
    code_dir = fileparts(code_file);
end
cd(code_dir);

figure_dir = fullfile(code_dir,'figures');
if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

K_ind = 4; Non_dimX = (1:2); Non_dimY = (1:2);
kappa_choice = 1; kappa_choice2 = 2;

data_file = fullfile(code_dir,'Data_FCRGT.xlsx');
temperature_data = readmatrix(data_file,'Sheet','LTemp');
grp_data = readmatrix(data_file,'Sheet','WGRP2');

GTemp_Grid = temperature_data(1,2:end);
GTemp0 = temperature_data(2:end,2:end);
T = size(GTemp0,1);

GRP_Grid = grp_data(1,2:end)';
GRP_Gr0 = grp_data(2:end,2:end);
Date = grp_data(2:end,1)';

%% Main Paper - Figure 1: Density Data and Descriptive Moments
fig_main1a = figure('Position',[100 100 1600 800],'Color','w');
mesh(GTemp_Grid',(1951:2019),GTemp0)
xlabel('Land Temperature Anomaly','Fontsize',15,'fontweight','b');
ylabel('Year','Fontsize',35,'fontweight','b')
xlim([GTemp_Grid(1)-0.05 GTemp_Grid(end)+0.05])
ylim([1951 2019])
set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main1a,fullfile(figure_dir,'GTemp_dist.png'),'Resolution',300,'BackgroundColor','white');

fig_main1b = figure('Position',[100 100 1600 800],'Color','w');
mesh(GRP_Grid,(1951:2019),GRP_Gr0)
xlabel('GRP Growth Rate','Fontsize',15,'fontweight','b');
ylabel('Year','Fontsize',35,'fontweight','b')
xlim([GRP_Grid(1)-0.05 GRP_Grid(end)+0.05])
ylim([1951 2019])
set(gca,'Linewidth',3.0,'box','off','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main1b,fullfile(figure_dir,'GRP_Gr_Dist_FE2.png'),'Resolution',300,'BackgroundColor','white');

Sum_stat_temp = KS_Desc_stat(GTemp_Grid',GTemp0);
fig_main1c = figure('Position',[100 100 1600 800],'Color','w');
plot(Date,Sum_stat_temp(:,1),'r','LineWidth',2.5); hold on;
plot(Date,sqrt(Sum_stat_temp(:,2)),'b','LineWidth',2.5); hold on;
axis tight; grid on;
legend('Mean','Standard Deviation','Location','Northwest','Orientation','horizontal'); legend boxoff;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main1c,fullfile(figure_dir,'GTemp_MV.png'),'Resolution',300,'BackgroundColor','white');

Sum_stat_grp = KS_Desc_stat(GRP_Grid,GRP_Gr0);
fig_main1d = figure('Position',[100 100 1600 800],'Color','w');
plot(Date,Sum_stat_grp(:,1),'r','LineWidth',2.5); hold on;
plot(Date,sqrt(Sum_stat_grp(:,2)),'b','LineWidth',2.5); hold on;
axis tight; grid on;
legend('Mean','Standard Deviation','Location','Northwest','Orientation','horizontal'); legend boxoff;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main1d,fullfile(figure_dir,'GRP_Gr_FE2.png'),'Resolution',300,'BackgroundColor','white');

%% Main Paper - Figure 2: Global Warming Functions GW1 and GW2
Xraw_tem = GTemp0'; T1 = size(Xraw_tem,2);
indx = floor(T1/2); GW2_size = 5;

fig_main2a = figure('Position',[100 100 1600 800],'Color','w');
plot(GTemp_Grid,mean(Xraw_tem(:,1:indx),2),'k--','LineWidth',3.5); hold on; grid on;
plot(GTemp_Grid,mean(Xraw_tem(:,(indx+1):T1),2),'k','LineWidth',3.5); hold on;
axis([GTemp_Grid(1) GTemp_Grid(end) -0.01 0.45]);
legend('1951-1984 Avg','1985-2019 Avg','Location','Northwest','box','off');
title('Global Warming Function (GW1)','fontsize',35,'fontweight','b','LineWidth',10.0);
xlabel('Temperature Anomaly','fontsize',35,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main2a,fullfile(figure_dir,'GW1_Shock.png'),'Resolution',300,'BackgroundColor','white');

fig_main2b = figure('Position',[100 100 1600 800],'Color','w');
plot(GTemp_Grid,mean(Xraw_tem(:,1:GW2_size),2),'k--','LineWidth',3.5); hold on; grid on;
plot(GTemp_Grid,mean(Xraw_tem(:,(T1-GW2_size+1):T1),2),'k','LineWidth',3.5); hold on;
axis([GTemp_Grid(1) GTemp_Grid(end) -0.01 0.45]);
legend('1951-1955 Avg','2015-2019 Avg','Location','Northwest','box','off');
title('Global Warming Function (GW2)','fontsize',35,'fontweight','b','LineWidth',10.0);
xlabel('Temperature Anomaly','fontsize',35,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main2b,fullfile(figure_dir,'GW2_Shock.png'),'Resolution',300,'BackgroundColor','white');

%% Shared Estimation for Main Figures 3-4 and Supplement Figures 1-5
GTemp = zeros(size(GTemp0)); GRP_Gr = zeros(size(GRP_Gr0));
for i = 1:T
    GTemp(i,:) = log(GTemp0(i,:)/geomean(GTemp0(i,:)));
    GRP_Gr(i,:) = log(GRP_Gr0(i,:)/geomean(GRP_Gr0(i,:)));
end

[~,Inv_N_GTemp,~,~,~,Inv_S_GTemp] = density_decom(GTemp_Grid,GTemp,max(Non_dimX));

E_GTemp = mean(GTemp)';
GTemp = GTemp-kron(ones(T,1),E_GTemp');
E_GRP_Gr = mean(GRP_Gr)';
GRP_Gr = GRP_Gr-kron(ones(T,1),E_GRP_Gr');

% Total-run estimators
Tef0_gw0 = f_kappa_est_ci(0,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);
Tef0_gw1 = f_kappa_est_ci(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);
Tef0_gw2 = f_kappa_est_ci(kappa_choice2,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);
Tef0_gw0_GW2 = f_kappa_est_ci_GW(0,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);
Tef0_gw1_GW2 = f_kappa_est_ci_GW(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);
Tef0_gw2_GW2 = f_kappa_est_ci_GW(kappa_choice2,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,1);

% Short-run estimators
Sef0_gw0 = f_kappa_est_ci(0,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);
Sef0_gw1 = f_kappa_est_ci(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);
Sef0_gw2 = f_kappa_est_ci(kappa_choice2,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);
Sef0_gw0_GW2 = f_kappa_est_ci_GW(0,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);
Sef0_gw1_GW2 = f_kappa_est_ci_GW(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);
Sef0_gw2_GW2 = f_kappa_est_ci_GW(kappa_choice2,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind,Non_dimY,Non_dimX,0);

%% Main Paper - Figure 3: Total-run and Short-run Response Functions
index = 2:(length(GRP_Grid)-1);

fig_main3a = figure('Position',[100 100 1600 800],'Color','w');
subplot(1,2,1)
plot(GRP_Grid(index),Tef0_gw0(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Tef0_gw1(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Tef0_gw1(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Tef0_gw1(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.45 0.45]); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Total-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0);
legend('GW1 at \kappa=0','GW1 at \kappa=1','95% CI','','NumColumns',1,'Location','Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw0(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Sef0_gw1(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Sef0_gw1(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.05 0.05]); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Short-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0);
legend('GW1 at \kappa=0','GW1 at \kappa=1','95% CI','','NumColumns',1,'Location','Southeast','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
exportgraphics(fig_main3a,fullfile(figure_dir,'SR_TR_Response_GW1.png'),'Resolution',300,'BackgroundColor','white');

fig_main3b = figure('Position',[100 100 1600 800],'Color','w');
subplot(1,2,1)
plot(GRP_Grid(index),Tef0_gw0_GW2(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Tef0_gw1_GW2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Tef0_gw1_GW2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Tef0_gw1_GW2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.8 0.8]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Total-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
legend('GW2 at \kappa=0','GW2 at \kappa=1','95% CI','','NumColumns',1,'Location','Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw0_GW2(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.15 0.2]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Short-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
legend('GW2 at \kappa=0','GW2 at \kappa=1','95% CI','','NumColumns',1,'Location','Southwest','box','off');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
exportgraphics(fig_main3b,fullfile(figure_dir,'SR_TR_Response_GW2.png'),'Resolution',300,'BackgroundColor','white');

%% Main Paper - Figure 4: Climate-impact Density Shifts and Moments
Wght = linspace(0,1.5,30);
Wght(7) = 0.3; Wght(13) = 0.6; Wght(19) = 0.9;
Wght(24) = 1.2; Wght(30) = 1.5;

[Target0_grid0,ref_recover,Proj_raw1,Proj_raw2,Proj_raw3,Proj_raw4,Proj_raw5,Sum_stat] = ...
    CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1,Wght);
[Target0_grid0_GW2,ref_recover_GW2,Proj_raw1_GW2,Proj_raw2_GW2,Proj_raw3_GW2,Proj_raw4_GW2,Proj_raw5_GW2,Sum_stat_GW2] = ...
    CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1_GW2,Wght);

fig_main4a = figure('Position',[100 100 1600 800],'Color','w');
plot(Target0_grid0(2:end-1),ref_recover(2:end-1),'k-.','LineWidth',4.5); hold on;
plot(Target0_grid0(2:end-1),Proj_raw1(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0(2:end-1),Proj_raw2(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0(2:end-1),Proj_raw3(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0(2:end-1),Proj_raw4(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0(2:end-1),Proj_raw5(2:end-1),'k','LineWidth',4.5); hold on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
axis([min(Target0_grid0) max(Target0_grid0) 0.2 1.7]); grid on;
legend('Avg GRP(51-84)','q=0.3','q=0.6','q=0.9','q=1.2','q=1.5','Location','Northoutside','box','off','Orientation','horizontal','Fontsize',40,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main4a,fullfile(figure_dir,'CC_Impact1.png'),'Resolution',300,'BackgroundColor','white');

fig_main4b = figure('Position',[100 100 1600 800],'Color','w');
plot(Target0_grid0_GW2(2:end-1),ref_recover_GW2(2:end-1),'k-.','LineWidth',4.5); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw1_GW2(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw2_GW2(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw3_GW2(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw4_GW2(2:end-1),'k:','LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw5_GW2(2:end-1),'k','LineWidth',4.5); hold on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
axis([min(Target0_grid0_GW2) max(Target0_grid0_GW2) 0.2 1.7]); grid on;
legend('Avg GRP(51-84)','q=0.3','q=0.6','q=0.9','q=1.2','q=1.5','Location','Northoutside','box','off','Orientation','horizontal','Fontsize',40,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_main4b,fullfile(figure_dir,'CC_Impact2.png'),'Resolution',300,'BackgroundColor','white');

fig_main4c = figure('Position',[100 100 1600 800],'Color','w');
subplot(2,1,1)
plot(Wght,Sum_stat(2:end,1),'b','LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2(2:end,1),'r','LineWidth',2.5); hold on;
title('Mean','fontsize',35,'fontweight','b','LineWidth',10.0); axis tight; grid on;
legend('GW1','GW2','Location','Northeast','box','off','Orientation','horizontal');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');
subplot(2,1,2)
plot(Wght,Sum_stat(2:end,2),'b','LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2(2:end,2),'r','LineWidth',2.5); hold on;
legend('GW1','GW2','Location','Northwest','box','off','Orientation','horizontal');
title('Variance','fontsize',35,'fontweight','b','LineWidth',10.0); axis tight; grid on;
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');
exportgraphics(fig_main4c,fullfile(figure_dir,'CC_Stat00.png'),'Resolution',300,'BackgroundColor','white');

%% Supplementary Material - Figure 1: Functional Residuals
[Resid1,Resid2,eval_S] = submission_residual_components( ...
    kappa_choice,GTemp_Grid,GTemp,GRP_Grid,GRP_Gr,K_ind,Non_dimX);

Num_basis = 200; Y_ngrid = length(GRP_Grid);
t_Y = (0:(Y_ngrid-1))/(Y_ngrid-1);
Func_Y = NaN(Y_ngrid,Num_basis);
for i = 1:(Num_basis/2)
    sin_func = sqrt(2)*sin(2*pi*i*t_Y);
    cos_func = sqrt(2)*cos(2*pi*i*t_Y);
    norm_sin = sqrt(inner_product(sin_func,sin_func,t_Y));
    norm_cos = sqrt(inner_product(cos_func,cos_func,t_Y));
    Func_Y(:,2*i-1) = sin_func/norm_sin;
    Func_Y(:,2*i) = cos_func/norm_cos;
end

Resid1_Func = (Func_Y*Resid1)'; Resid2_Func = (Func_Y*Resid2)';
Resid_Date = Date((kappa_choice+1):end);
fig_supp1 = figure('Position',[100 100 2000 1000],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
[X_resid,Y_resid] = meshgrid(GRP_Grid,Resid_Date);
clim_val = [-2 2]; ncol = 256;
cmap = interp1([1,ncol/2,ncol],[0.00 0.10 0.65;0.45 0.70 0.95;0.70 0.00 0.00],1:ncol);

nexttile;
mesh(X_resid,Y_resid,Resid1_Func,Resid1_Func,'EdgeColor','interp','LineWidth',2.5);
colormap(cmap); clim(clim_val); zlim([-2 2]); view(45,25);
xlabel('GRP Gr','FontSize',33,'FontWeight','bold');
ylabel('Year','FontSize',33,'FontWeight','bold');
zlabel('Residual','FontSize',33,'FontWeight','bold');
title('$y_t-\widehat{f}^{N}_{\kappa}(\widetilde{x}_t)$','Interpreter','latex','FontSize',48,'FontWeight','bold');
xlim([GRP_Grid(1)-0.005 GRP_Grid(end)+0.005]); ylim([Resid_Date(1) Resid_Date(end)]);
set(gca,'LineWidth',3.5,'Box','off','TickLength',[0.01 0.01],'FontSize',33,'FontWeight','bold','FontName','Arial');
pbaspect([1 1.6 0.8]); cb1 = colorbar; cb1.FontSize = 33; cb1.FontWeight = 'bold';

nexttile;
mesh(X_resid,Y_resid,Resid2_Func,Resid2_Func,'EdgeColor','interp','LineWidth',2.5);
colormap(cmap); clim(clim_val); zlim([-2 2]); view(45,25);
xlabel('GRP Gr','FontSize',33,'FontWeight','bold');
ylabel('Year','FontSize',33,'FontWeight','bold');
zlabel('Residual','FontSize',33,'FontWeight','bold');
title('$y_t-\widehat{f}_{\kappa}(\widetilde{x}_t)$','Interpreter','latex','FontSize',48,'FontWeight','bold');
xlim([GRP_Grid(1)-0.005 GRP_Grid(end)+0.005]); ylim([Resid_Date(1) Resid_Date(end)]);
set(gca,'LineWidth',3.5,'Box','off','TickLength',[0.01 0.01],'FontSize',33,'FontWeight','bold','FontName','Arial');
pbaspect([1 1.6 0.8]); cb2 = colorbar; cb2.FontSize = 33; cb2.FontWeight = 'bold';
exportgraphics(fig_supp1,fullfile(figure_dir,'FCRGT_Residuals.png'),'Resolution',600,'BackgroundColor','white');

%% Supplementary Material - Figure 2: Stationary-dimension Diagnostics
scree_index = 1:min(10,length(eval_S));

fig_supp2a = figure('Position',[100 100 900 650],'Color','w');
plot(scree_index,eval_S(scree_index),'k-o','LineWidth',3.0,'MarkerSize',8); hold on;
grid on; axis tight;
xlabel('Stationary Component','fontsize',20,'fontweight','b');
ylabel('Eigenvalue','fontsize',20,'fontweight','b');
title('Stationary-part Scree Plot','fontsize',22,'fontweight','b','LineWidth',10.0);
legend('Eigenvalue','Location','Northeast','box','off','FontSize',16);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'XTick',scree_index,'Fontsize',18,'fontweight','b');
exportgraphics(fig_supp2a,fullfile(figure_dir,'FCRGT_Stationary_Scree.png'),'Resolution',300,'BackgroundColor','white');

eval_S_share = 100*eval_S/sum(eval_S); eval_S_cumshare = cumsum(eval_S_share);
fig_supp2b = figure('Position',[100 100 900 650],'Color','w');
b_share = bar(scree_index,eval_S_share(scree_index),0.65,'FaceColor',[0.35 0.65 0.90],'EdgeColor','none'); hold on;
p_cum = plot(scree_index,eval_S_cumshare(scree_index),'k-o','LineWidth',3.0,'MarkerSize',8);
h95 = yline(95,'r--','LineWidth',2.5); grid on;
text(1,eval_S_share(1)+3,[num2str(eval_S_share(1),'%.1f') '%'],'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
text(2,eval_S_share(2)+3,[num2str(eval_S_share(2),'%.1f') '%'],'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
xlabel('Stationary Component','fontsize',20,'fontweight','b');
ylabel('Explained Proportion (%)','fontsize',20,'fontweight','b');
title('Explained Proportions of Stationary Components','fontsize',22,'fontweight','b','LineWidth',10.0);
legend([b_share p_cum h95],{'Individual','Cumulative','95% threshold'},'Location','Southeast','box','off','FontSize',15);
axis([0.5 max(scree_index)+0.5 0 105]);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'XTick',scree_index,'YTick',0:10:100,'Fontsize',18,'fontweight','b');
exportgraphics(fig_supp2b,fullfile(figure_dir,'FCRGT_Stationary_Explained_Share.png'),'Resolution',300,'BackgroundColor','white');

%% Supplementary Material - Figure 3: Robustness to Kappa = 2
fig_supp3a = figure('Position',[100 100 1600 650],'Color','w');
subplot(1,2,1)
plot(GRP_Grid(index),Tef0_gw1(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Tef0_gw2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Tef0_gw2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Tef0_gw2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.45 0.45]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Total-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0);
legend('GW1 at \kappa=1','GW1 at \kappa=2','95% CI','','NumColumns',1,'Location','Southwest','box','on','Color','white','EdgeColor','white','FontSize',35);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw1(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Sef0_gw2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Sef0_gw2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.05 0.05]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Short-run Response Function','fontsize',22,'fontweight','b','LineWidth',10.0);
legend('GW1 at \kappa=1','GW1 at \kappa=2','95% CI','','NumColumns',1,'Location','Southeast','box','on','Color','white','EdgeColor','white','FontSize',35);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
exportgraphics(fig_supp3a,fullfile(figure_dir,'FCRGT_Kappa2_GW1.png'),'Resolution',300,'BackgroundColor','white');

fig_supp3b = figure('Position',[100 100 1600 650],'Color','w');
subplot(1,2,1)
plot(GRP_Grid(index),Tef0_gw1_GW2(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Tef0_gw2_GW2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Tef0_gw2_GW2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Tef0_gw2_GW2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.8 0.8]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Total-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0);
legend('GW2 at \kappa=1','GW2 at \kappa=2','95% CI','','NumColumns',1,'Location','Southwest','box','on','Color','white','EdgeColor','white','FontSize',35);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw1_GW2(index,1),'r','LineWidth',2.0); hold on;
plot(GRP_Grid(index),Sef0_gw2_GW2(index,1),'k','LineWidth',4.0); hold on;
plot(GRP_Grid(index),Sef0_gw2_GW2(index,2),'k--','LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw2_GW2(index,3),'k--','LineWidth',3.0); hold on;
axis([min(GRP_Grid) max(GRP_Grid) -0.15 0.2]); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Short-run Response Function','fontsize',35,'fontweight','b','LineWidth',10.0);
legend('GW2 at \kappa=1','GW2 at \kappa=2','95% CI','','NumColumns',1,'Location','Southwest','box','on','Color','white','EdgeColor','white','FontSize',35);
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
exportgraphics(fig_supp3b,fullfile(figure_dir,'FCRGT_Kappa2_GW2.png'),'Resolution',300,'BackgroundColor','white');

%% Supplementary Material - Figure 4: Robustness to K = 3, 5, and 6
K_ind_robust = [3 5 6];
Tef0_gw1_robust = zeros(length(GRP_Grid),3,length(K_ind_robust));
Tef0_gw1_GW2_robust = zeros(length(GRP_Grid),3,length(K_ind_robust));
Sef0_gw1_robust = zeros(length(GRP_Grid),3,length(K_ind_robust));
Sef0_gw1_GW2_robust = zeros(length(GRP_Grid),3,length(K_ind_robust));
for iK = 1:length(K_ind_robust)
    K_ind_alt = K_ind_robust(iK);
    Tef0_gw1_robust(:,:,iK) = f_kappa_est_ci(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind_alt,Non_dimY,Non_dimX,1);
    Tef0_gw1_GW2_robust(:,:,iK) = f_kappa_est_ci_GW(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind_alt,Non_dimY,Non_dimX,1);
    Sef0_gw1_robust(:,:,iK) = f_kappa_est_ci(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind_alt,Non_dimY,Non_dimX,0);
    Sef0_gw1_GW2_robust(:,:,iK) = f_kappa_est_ci_GW(kappa_choice,GTemp_Grid,GTemp0,GTemp,GRP_Grid,GRP_Gr0,GRP_Gr,Inv_N_GTemp,Inv_S_GTemp,K_ind_alt,Non_dimY,Non_dimX,0);
end
K_color = [0.0000 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880];

fig_supp4a = figure('Position',[100 100 1600 650],'Color','w');
subplot(1,2,1)
h_ci = plot(GRP_Grid(index),Tef0_gw1(index,2),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Tef0_gw1(index,3),'k--','LineWidth',2.5); hold on;
h_k3 = plot(GRP_Grid(index),Tef0_gw1_robust(index,1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
h_k5 = plot(GRP_Grid(index),Tef0_gw1_robust(index,1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
h_k6 = plot(GRP_Grid(index),Tef0_gw1_robust(index,1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
h_k4 = plot(GRP_Grid(index),Tef0_gw1(index,1),'k','LineWidth',4.0); hold on;
xlim([min(GRP_Grid) max(GRP_Grid)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Total-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
legend([h_k4 h_k3 h_k5 h_k6 h_ci],{'K=4 (baseline)','K=3','K=5','K=6','95% CI for K=4'},'NumColumns',1,'Location','Northeast','box','off','FontSize',18);
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw1(index,2),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Sef0_gw1(index,3),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Sef0_gw1_robust(index,1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_robust(index,1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_robust(index,1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1(index,1),'k','LineWidth',4.0); hold on;
xlim([min(GRP_Grid) max(GRP_Grid)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
title('Short-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
exportgraphics(fig_supp4a,fullfile(figure_dir,'FCRGT_K_Robustness_GW1.png'),'Resolution',300,'BackgroundColor','white');

fig_supp4b = figure('Position',[100 100 1600 650],'Color','w');
subplot(1,2,1)
h_ci = plot(GRP_Grid(index),Tef0_gw1_GW2(index,2),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Tef0_gw1_GW2(index,3),'k--','LineWidth',2.5); hold on;
h_k3 = plot(GRP_Grid(index),Tef0_gw1_GW2_robust(index,1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
h_k5 = plot(GRP_Grid(index),Tef0_gw1_GW2_robust(index,1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
h_k6 = plot(GRP_Grid(index),Tef0_gw1_GW2_robust(index,1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
h_k4 = plot(GRP_Grid(index),Tef0_gw1_GW2(index,1),'k','LineWidth',4.0); hold on;
xlim([min(GRP_Grid) max(GRP_Grid)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Total-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
legend([h_k4 h_k3 h_k5 h_k6 h_ci],{'K=4 (baseline)','K=3','K=5','K=6','95% CI for K=4'},'NumColumns',1,'Location','Northeast','box','off','FontSize',18);
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
subplot(1,2,2)
plot(GRP_Grid(index),Sef0_gw1_GW2(index,2),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2(index,3),'k--','LineWidth',2.5); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2_robust(index,1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2_robust(index,1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2_robust(index,1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
plot(GRP_Grid(index),Sef0_gw1_GW2(index,1),'k','LineWidth',4.0); hold on;
xlim([min(GRP_Grid) max(GRP_Grid)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',30,'fontweight','b');
title('Short-run Response Function','fontsize',30,'fontweight','b','LineWidth',10.0);
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
exportgraphics(fig_supp4b,fullfile(figure_dir,'FCRGT_K_Robustness_GW2.png'),'Resolution',300,'BackgroundColor','white');

%% Supplementary Material - Figure 5: Climate-impact Robustness to K
Proj_raw5_GW1_robust = zeros(length(GRP_Grid),length(K_ind_robust));
Proj_raw5_GW2_robust = zeros(length(GRP_Grid),length(K_ind_robust));
Sum_stat_GW1_robust = zeros(length(Wght)+1,4,length(K_ind_robust));
Sum_stat_GW2_robust = zeros(length(Wght)+1,4,length(K_ind_robust));
for iK = 1:length(K_ind_robust)
    [~,~,~,~,~,~,Proj_raw5_alt,Sum_stat_alt] = CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1_robust(:,:,iK),Wght);
    Proj_raw5_GW1_robust(:,iK) = Proj_raw5_alt;
    Sum_stat_GW1_robust(:,:,iK) = Sum_stat_alt;
    [~,~,~,~,~,~,Proj_raw5_alt,Sum_stat_alt] = CC_Frac_Impacts(GRP_Grid,GRP_Gr0,Tef0_gw1_GW2_robust(:,:,iK),Wght);
    Proj_raw5_GW2_robust(:,iK) = Proj_raw5_alt;
    Sum_stat_GW2_robust(:,:,iK) = Sum_stat_alt;
end

fig_supp5a = figure('Position',[100 100 2400 900],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax_CC_K_density1 = nexttile;
h_ref = plot(Target0_grid0(2:end-1),ref_recover(2:end-1),'k-.','LineWidth',4.5); hold on;
h_k3 = plot(Target0_grid0(2:end-1),Proj_raw5_GW1_robust(2:end-1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
h_k5 = plot(Target0_grid0(2:end-1),Proj_raw5_GW1_robust(2:end-1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
h_k6 = plot(Target0_grid0(2:end-1),Proj_raw5_GW1_robust(2:end-1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
h_k4 = plot(Target0_grid0(2:end-1),Proj_raw5(2:end-1),'k','LineWidth',4.5); hold on;
xlim([min(Target0_grid0) max(Target0_grid0)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b'); title('GW1 at q=1.5','fontsize',35,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');

nexttile;
plot(Target0_grid0_GW2(2:end-1),ref_recover_GW2(2:end-1),'k-.','LineWidth',4.5); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw5_GW2_robust(2:end-1,1),'Color',K_color(1,:),'LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw5_GW2_robust(2:end-1,2),'Color',K_color(2,:),'LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw5_GW2_robust(2:end-1,3),'Color',K_color(3,:),'LineWidth',3.0); hold on;
plot(Target0_grid0_GW2(2:end-1),Proj_raw5_GW2(2:end-1),'k','LineWidth',4.5); hold on;
xlim([min(Target0_grid0_GW2) max(Target0_grid0_GW2)]); ylim('padded'); grid on;
xlabel('GRP Growth Rate','fontsize',35,'fontweight','b'); title('GW2 at q=1.5','fontsize',35,'fontweight','b');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
lgd_CC_K_density = legend(ax_CC_K_density1,[h_ref h_k4 h_k3 h_k5 h_k6],{'Avg GRP (51-84)','K=4','K=3','K=5','K=6'},'Orientation','horizontal','NumColumns',5,'box','off','FontSize',40,'FontWeight','bold');
lgd_CC_K_density.Layout.Tile = 'north';
exportgraphics(fig_supp5a,fullfile(figure_dir,'FCRGT_K_Robustness_Climate_Impact_Densities.png'),'Resolution',300,'BackgroundColor','white');

fig_supp5b = figure('Position',[100 100 1500 1200],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
h_k4 = plot(Wght,Sum_stat(2:end,1),'k','LineWidth',3.0); hold on;
plot(Wght,Sum_stat_GW2(2:end,1),'k--','LineWidth',3.0); hold on;
h_k3 = plot(Wght,Sum_stat_GW1_robust(2:end,1,1),'Color',K_color(1,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,1,1),'--','Color',K_color(1,:),'LineWidth',2.5); hold on;
h_k5 = plot(Wght,Sum_stat_GW1_robust(2:end,1,2),'Color',K_color(2,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,1,2),'--','Color',K_color(2,:),'LineWidth',2.5); hold on;
h_k6 = plot(Wght,Sum_stat_GW1_robust(2:end,1,3),'Color',K_color(3,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,1,3),'--','Color',K_color(3,:),'LineWidth',2.5); hold on;
xlim([min(Wght) max(Wght)]); ylim('padded'); grid on; title('Mean','fontsize',35,'fontweight','b');
legend([h_k4 h_k3 h_k5 h_k6],{'K=4','K=3','K=5','K=6'},'Location','Southwest','NumColumns',2,'box','off','FontSize',33,'FontWeight','bold');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');

nexttile;
h_gw1 = plot(Wght,Sum_stat(2:end,2),'k','LineWidth',3.0); hold on;
h_gw2 = plot(Wght,Sum_stat_GW2(2:end,2),'k--','LineWidth',3.0); hold on;
plot(Wght,Sum_stat_GW1_robust(2:end,2,1),'Color',K_color(1,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,2,1),'--','Color',K_color(1,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW1_robust(2:end,2,2),'Color',K_color(2,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,2,2),'--','Color',K_color(2,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW1_robust(2:end,2,3),'Color',K_color(3,:),'LineWidth',2.5); hold on;
plot(Wght,Sum_stat_GW2_robust(2:end,2,3),'--','Color',K_color(3,:),'LineWidth',2.5); hold on;
xlim([min(Wght) max(Wght)]); ylim('padded'); grid on; title('Variance','fontsize',35,'fontweight','b');
legend([h_gw1 h_gw2],{'GW1 (solid)','GW2 (dashed)'},'Location','Northwest','box','off','FontSize',33,'FontWeight','bold');
set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',33,'fontweight','b');
exportgraphics(fig_supp5b,fullfile(figure_dir,'FCRGT_K_Robustness_Climate_Impact_Moments.png'),'Resolution',300,'BackgroundColor','white');

fprintf('All paper-linked figures were written to:\n%s\n',figure_dir);

%% Local Function: Paper-relevant Residual Components Only
function [Resid1,Resid2,eval_S] = submission_residual_components( ...
    kappa,GTemp_Grid,GTemp,GRP_Grid,GRP_Gr,K_ind,Non_dimX)

T = size(GTemp,1); Num_basis = 200;

X_ngrid = length(GTemp_Grid);
t_X = (0:(X_ngrid-1))/(X_ngrid-1);
Func_X = NaN(X_ngrid,Num_basis);
for i = 1:(Num_basis/2)
    sin_func = sqrt(2)*sin(2*pi*i*t_X);
    cos_func = sqrt(2)*cos(2*pi*i*t_X);
    norm_sin = sqrt(inner_product(sin_func,sin_func,t_X));
    norm_cos = sqrt(inner_product(cos_func,cos_func,t_X));
    Func_X(:,2*i-1) = sin_func/norm_sin;
    Func_X(:,2*i) = cos_func/norm_cos;
end
XX = (Func_X'*GTemp')*(t_X(2)-t_X(1));

Y_ngrid = length(GRP_Grid);
t_Y = (0:(Y_ngrid-1))/(Y_ngrid-1);
Func_Y = NaN(Y_ngrid,Num_basis);
for i = 1:(Num_basis/2)
    sin_func = sqrt(2)*sin(2*pi*i*t_Y);
    cos_func = sqrt(2)*cos(2*pi*i*t_Y);
    norm_sin = sqrt(inner_product(sin_func,sin_func,t_Y));
    norm_cos = sqrt(inner_product(cos_func,cos_func,t_Y));
    Func_Y(:,2*i-1) = sin_func/norm_sin;
    Func_Y(:,2*i) = cos_func/norm_cos;
end
YY = (Func_Y'*GRP_Gr')*(t_Y(2)-t_Y(1));

FC_X1 = XX(:,1:(T-kappa));
FC_X0 = XX(:,(kappa+1):T);
FC_Y = YY(:,(kappa+1):T);

C_kap = (FC_X0*FC_X1')/T;
D_kap = C_kap'*C_kap;
[Y_eigvecs,Y_eigvals] = eig(D_kap);
[eval_D,Y_indsrt] = sort(diag(Y_eigvals),'descend');
evec_D = Y_eigvecs(:,Y_indsrt);

Z1 = (FC_X1'*evec_D(:,1:K_ind))';
Z0 = (FC_X0'*evec_D(:,1:K_ind))';
D_kap_inv = diag(1./eval_D(1:K_ind));
C_kap_D = (Z0*Z1')/T;
CR_kap = (FC_Y*Z1')/T;

diag_N = [ones(max(Non_dimX),1);zeros(K_ind-max(Non_dimX),1)];
fkap_N = CR_kap*C_kap_D*D_kap_inv*diag(diag_N);
FC_resid = FC_Y-fkap_N*Z0;

CR_kap2 = (FC_resid*Z1')/T;
diag_S = [zeros(max(Non_dimX),1);ones(K_ind-max(Non_dimX),1)];
fkap_S = CR_kap2*C_kap_D*D_kap_inv*diag(diag_S);

Resid1 = FC_resid;
Resid2 = FC_Y-(fkap_N+fkap_S)*Z0;

Cov_S = (FC_resid*FC_resid')/size(FC_resid,2);
Cov_S = (Cov_S+Cov_S')/2;
eval_S = sort(real(eig(Cov_S)),'descend');
end
