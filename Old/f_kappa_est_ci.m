function [Tef_gw_store] = f_kappa_est_ci(kappa, GTemp_Grid, GTemp0, GTemp, GRP_Grid, GRP_Gr0, GRP_Gr, Inv_N_GTemp, Inv_S_GTemp, K_ind, Non_dimY, Non_dimX, TR)
%#####################################################
%######## End of 1. PREPARATION ######################
%#####################################################  
T = size(GTemp,1);  Num_basis=200;
X_ngrid = length(GTemp_Grid); t_X = (0:(X_ngrid-1))/(X_ngrid-1);

% Initialize the basis function matrix LBF
Func_X = NaN(X_ngrid, Num_basis);

% Compute the sine and cosine basis functions
for i = 1:(Num_basis/2)
    sin_func = sqrt(2)*sin(2*pi*i*t_X); cos_func = sqrt(2)*cos(2*pi*i*t_X);
    norm_sin = sqrt(inner_product(sin_func,sin_func,t_X)); norm_cos = sqrt(inner_product(cos_func,cos_func,t_X));
    Func_X(:,2*i-1) = sin_func/norm_sin; Func_X(:,2*i) = cos_func/norm_cos;
end

Xraw = GTemp; % Assuming GTemp is the input data
X = Xraw'; % Transpose for matrix operations

% Compute XX using the basis functions and the normalized time step
XX = (Func_X(1:X_ngrid,:)'*X(1:X_ngrid,:))*(t_X(2)-t_X(1));
%####################################################
% Representation of Y using Non_dimY
%####################################################
Y_ngrid = length(GRP_Grid); t_Y = (0:(Y_ngrid-1))/(Y_ngrid-1);

% Initialize the basis function matrix for Y
Func_Y = NaN(Y_ngrid, Num_basis);

% Compute the sine and cosine basis functions for Y
for i = 1:(Num_basis / 2)
    sin_func = sqrt(2)*sin(2*pi*i*t_Y); cos_func = sqrt(2)*cos(2*pi*i*t_Y);
    norm_sin = sqrt(inner_product(sin_func,sin_func,t_Y)); norm_cos = sqrt(inner_product(cos_func,cos_func,t_Y));
    Func_Y(:, 2*i-1) = sin_func/norm_sin; Func_Y(:,2*i) = cos_func/norm_cos;
end

% Process Yraw (GRP_Gr) with the basis functions
Yraw = GRP_Gr; % Assuming GRP_Gr is the input data
Y = Yraw'; % Transpose for matrix operations

% Compute YY using the basis functions and the normalized time step
YY = (Func_Y(1:Y_ngrid,:)'*Y(1:Y_ngrid,:))*(t_Y(2)-t_Y(1));

%#####################################################
%######## 2. ESTIMATION ##############################
%#####################################################
% Functional Representation with Trigonometric Basis are contained in XX and YY
FC_X1 = XX(:,1:(T-kappa));
FC_X0 = XX(:,(kappa+1):T);
FC_Y = YY(:,(kappa+1):T);

% Eigen decomposition of covariance matrices
C_kap = (FC_X0*FC_X1')/T;
D_kap = C_kap'*C_kap;

[Y_eigvecs, Y_eigvals] = eig(D_kap); % Eigenvalues and Eigenvectors
[eval_D, Y_indsrt] = sort(diag(Y_eigvals),'descend'); % Ordered eigenvalues
evec_D = Y_eigvecs(:,Y_indsrt); % Ordered eigenvectors

% Z matrices for estimation
Z1 = (FC_X1'*evec_D(:,1:K_ind))';
Z0 = (FC_X0'*evec_D(:,1:K_ind))';
Z1S = (FC_X1'*evec_D(:,(max(Non_dimX)+1):K_ind))';
Z0S = (FC_X0'*evec_D(:,(max(Non_dimX)+1):K_ind))';

% Adjusted covariance matrices in reduced space
C_kap_D = (Z0*Z1')/T;
D_kap_inv = diag(1./eval_D(1:K_ind));
CR_kap = (FC_Y*Z1')/T;

% f_N estimator
diag_N = [ones(max(Non_dimX), 1); zeros(K_ind-max(Non_dimX), 1)];
fkap_N = CR_kap*C_kap_D*D_kap_inv*diag(diag_N);

% Stationary residuals for f_S estimator
FC_resid = FC_Y - fkap_N*Z0;

CR_kap2 = (FC_resid * Z1')/T;

% f_S estimator
diag_S = [zeros(max(Non_dimX), 1); ones(K_ind - max(Non_dimX), 1)];
fkap_S = CR_kap2 * C_kap_D * D_kap_inv * diag(diag_S);

% Residuals from cointegrating regression
RESID = FC_Y - (fkap_N + fkap_S) * (evec_D(:, 1:K_ind)') * FC_X0;

%########################################################
%##### STEP2: Inference on Global Warming (gw) #########
%########################################################
%%% Generate the distributional Shock for Global Warming
Xraw_tem = GTemp0';
T1 = size(Xraw_tem, 2);
% indx = floor(T1/2)-7; M = 16; % 51-77 vs 93 - 2019
indx = floor(T1 / 2)-0; M = 1; % 51-78 vs 92 - 2019

gw_f = mean(Xraw_tem(:,(indx+M):T1),2) - mean(Xraw_tem(:,1:indx),2);
gw_clr = log(mean(Xraw_tem(:,(indx+M):T1),2) ./ geomean(mean(Xraw_tem(:,(indx+M):T1),2))) - ...
         log(mean(Xraw_tem(:,1:indx),2) ./ geomean(mean(Xraw_tem(:,1:indx),2)));
gw = Func_X(1:X_ngrid, :)'*gw_clr * (t_X(2) - t_X(1));

Nonstat_X_GW = zeros(T,1); Stat_X_GW = zeros(T,1); 
for i = 1:T
    Nonstat_X_GW(i) = mrsum(GTemp_Grid',(gw_f.*Inv_N_GTemp(i,:)')); 
    Stat_X_GW(i) = mrsum(GTemp_Grid',(gw_f.*Inv_S_GTemp(i,:)')); 
end

% figure;
% subplot(1,2,1); 
% plot((1951:1:2019),Nonstat_X_GW,'r','LineWidth',2.5); hold on; axis tight; grid on;
% xlabel('Year','fontsize',20,'fontweight','b'); legend('Nonstat- GW', 'Location', 'Northwest','box','off');
% set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b'); 
% subplot(1,2,2); 
% plot((1951:1:2019),Stat_X_GW,'r','LineWidth',2.5); hold on; axis tight; grid on;
% xlabel('Year','fontsize',20,'fontweight','b'); legend('Stat- GW', 'Location', 'Northwest','box','off');
% set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',30,'fontweight','b');
%%
if TR == 1;
e_gw = (fkap_N + fkap_S) * (evec_D(:, 1:K_ind)') * gw;
elseif TR == 0;
e_gw = (fkap_S) * (evec_D(:, 1:K_ind)') * gw;
end

Tef_gw = Func_Y*e_gw;
SD_kap_inv = diag(1 ./ eval_D((max(Non_dimX) + 1):K_ind));
SC_kap_D = (Z0S * Z1S') / T;
SC_D = (Z0S * Z0S') / T;

P_KS = evec_D(:,(max(Non_dimX)+1):K_ind)';
theta_raw = (SD_kap_inv * (SC_kap_D' * SC_D * SC_kap_D * SD_kap_inv) * P_KS * gw) .* (P_KS * gw);
% theta = sqrt(sum(theta_raw .^ 2));
theta = (sum(theta_raw));
C_resid = (RESID * RESID') / T;

stepsize = 1;
III = [];
SDFAC = [];
for i = 1:length(GRP_Grid)
    II = zeros(length(GRP_Grid), 1);
    II(max(i - stepsize, 1):min(i + stepsize, length(GRP_Grid))) = 1;
    III = [III; II' / sum(II)];
    IIF = Func_Y' * II;
    SDFAC = [SDFAC; 1.96 * sqrt(theta * (IIF' * C_resid * IIF) / T)];
end

Tef_gwu = Tef_gw + SDFAC;
Tef_gwl = Tef_gw - SDFAC;

lmu = max(Tef_gwu);
lml = min(Tef_gwl);

index = 2:(length(GRP_Grid)-1);
% figure;
% plot(GRP_Grid(index), Tef_gw(index), 'r','LineWidth',2.0); hold on;
% plot(GRP_Grid(index), Tef_gwu(index), 'k--','LineWidth',2.0); hold on;
% plot(GRP_Grid(index), Tef_gwl(index), 'k--','LineWidth',2.0); hold on;
% axis tight; grid on;
% xlabel('GRP Growth Rate','fontsize',35,'fontweight','b');
% title(['Global Warming with \kappa = ', num2str(kappa)],'fontsize',35,'fontweight','b','LineWidth',10.0); 
% legend('TR\_GW', '95% CI', '', 'Location', 'Northeast','box','off');
% set(gca,'Linewidth',3.0,'box','on','Ticklength',[0 0],'Fontsize',35,'fontweight','b');
Tef_gw_store = [Tef_gw Tef_gwu Tef_gwl];



