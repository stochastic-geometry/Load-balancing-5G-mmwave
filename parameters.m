% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Author: Chiranjib Saha, Harpreet S. Dhillon
% Email: csaha@vt.edu, hdhillon@vt.edu
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% This code was used to make coverage plots in the following paper
% Bibentry goes here ----
%- ------


% Specify all parameter values in this script


% This script is written to test rate coverage with bias factor sweep
%% Spatial distribution parametes
% Ne-6 where N is number of points per sq km
  lambda_UE = 1000e-6; % UE PPP intensity
% % % lambda_UE = 0 % set this to zero to speed up simulation for coverage probability (debug).
% % % 
 lambda_MBS = 5e-6;  % MBS PPP intensity
 if ~exist('lambda_SBS')
  lambda_SBS = 100e-6; % SBS PPP intensity
 end
%% Simulation parameters
diskRadius = 1000; % set so that there are on average 30 MBSs in simulation window
diskRadius_UE = diskRadius;
MaxIter  = 5e3;
diskArea=pi*diskRadius^2;
diskArea_UE=pi*diskRadius_UE^2;
NumPoints = 15;
if ~exist('Rate_threshold')
 Rate_threshold = 1.937e7;%50e6;%logspace(6.5,9,NumPoints);
end
 %% Pathloss 
alpha_l = 3.0;
alpha_n = 4.0;
L = @(x,state) (state== 1).*x.^-alpha_l + (state==0).*x.^-alpha_n;
lambda_carrier = 3e8/28e9;
beta = (lambda_carrier/(4*pi))^2;
betaConst=beta; % this is to avoid confict of predefined `beta' function in Matlab
%% Shadowing
Shadow_coefficient = 4; % dB
is_shadowing = false;
%% Transmit powers and Bias Factors 
P_m = 10^(40/10);
%P_m = 10^(24/20);
P_s = 10^(20/10);

T_s_db =0;%linspace(0,50,NumPoints);
T_s_lin = 10.^(T_s_db./10);   % 10;
T_s = T_s_lin;
T_m = 1;%/B_s;
if ~exist('cov_th')
 cov_th =1;
end
cov_th_b = 10^(5/10);
%% System-level parameters
W =  1000e6; % System BW
N0 = 10^((10*log10(W)-173+10)/10); % noise PSD
N_0 = N0/beta; % in analysis, beta is not used explicitly in the formulas. 

%% Set output filenames
outfile = 'Simualtion_part0';
%% Antenna Parameters
G_m = 10^(18/10);
G_s = G_m;
g_m = 10^(-2/10);
g_s = g_m; 
theta_hpbw= 10*pi/180;
%% Blockage
mu = 200;
%%%%%%%%%%%%%%%%%%%%%%

eta = 0.8; %bandwidth partition factor for ORA
%% Specifications for Germ Grain model

% MU = PI/(Lambda_Bl L_BL)

lambda_Block =  15e-04;
Block_length = 5;
