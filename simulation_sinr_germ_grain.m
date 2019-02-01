
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Author: Chiranjib Saha, Harpreet S. Dhillon
% Email: csaha@vt.edu, hdhillon@vt.edu
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% This code was used to make coverage plots in the following paper
% Bibentry goes here ----
%- ------


% Run this script to generate the association probability for a given
% configuration (lambda_bl,L_bl) of blockage distribution. Set these values
% in "parameters.m". 



clear all; close all;
parameters;

fprintf('\n Starting simulation...\n');
% initialize counter
internal_counter_m = 0;
for count_sim=1:MaxIter
  %% Generate number of MBSs, SBSs, and blockages in the simulation disk    
    randNumb_MBS=poissrnd(lambda_MBS*diskArea);
    randNumb_SBS=poissrnd(lambda_SBS*diskArea);
    randNumb_Block=poissrnd(lambda_Block*diskArea);
    
  %% step 1: generate MBS PPPs
  % Generating MBS PPP
    theta = rand(randNumb_MBS,1)*(2*pi);
    r = diskRadius*sqrt(rand(randNumb_MBS,1));
    x =  r.*cos(theta);  
    y =  r.*sin(theta);   
    MBS_location=[x,y];

  % Generating SBS PPP
    theta = rand(randNumb_SBS,1)*(2*pi);
    r = diskRadius*sqrt(rand(randNumb_SBS,1));
    x =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
    y =  r.*sin(theta);   %%%************************************************
    SBS_location=[x,y];
  % Generating Blockage PPP
    theta = rand(randNumb_Block,1)*(2*pi);
    r = diskRadius*sqrt(rand(randNumb_Block,1));
    x =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
    y =  r.*sin(theta);   %%%************************************************
    Block_location=[x,y];  
    Block_Orientation = rand(randNumb_Block,1)*(2*pi);
    Block_endpoints_x  = [Block_location(:,1)-Block_length/2*cos(Block_Orientation),Block_location(:,1)+Block_length/2*cos(Block_Orientation)];
    Block_endpoints_y  = [Block_location(:,2)-Block_length/2*sin(Block_Orientation),Block_location(:,2)+Block_length/2*sin(Block_Orientation)];

% Generating UE PPP

   UE_location=[0,0];  
   
  %%%% Comment this out while running the full simulation %% 
  %%%% Plot the network for visulization purpose %%%% 
%    plot([Block_endpoints_x'],[Block_endpoints_y'],'r-','linewidth',2)
%    axis('square');
%    hold on;
%    plot(MBS_location(:,1),MBS_location(:,2),'o','linewidth',2);
%    plot(UE_location(:,1),UE_location(:,2),'k.');
%    hold off;
%    keyboard;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%
  %% Step 2: Link state Computation 
  % First construct the LinkState matrix
  % For MBS-UE
    MBS_loc_rep = repmat(MBS_location,size(UE_location,1),1);
    UE_loc_rep  = repelem(UE_location,size(MBS_location,1),1);
    out= lineSegmentIntersect([UE_loc_rep,MBS_loc_rep],[Block_endpoints_x(:,1),Block_endpoints_y(:,1),...
                  Block_endpoints_x(:,2),Block_endpoints_y(:,2)]);
   Linkstate_MBS_UE_int = reshape(sum(out,2),size(MBS_location,1),size(UE_location,1));
   Linkstate_MBS_UE = (Linkstate_MBS_UE_int==0);
   Linkstate_MBS_UE = Linkstate_MBS_UE';  % for dimension matching in later steps
  % For SBS-UE
   SBS_loc_rep = repmat(SBS_location,size(UE_location,1),1);
   UE_loc_rep  = repelem(UE_location,size(SBS_location,1),1);
   out= lineSegmentIntersect([UE_loc_rep,SBS_loc_rep],[Block_endpoints_x(:,1),Block_endpoints_y(:,1),...
                  Block_endpoints_x(:,2),Block_endpoints_y(:,2)]);
   Linkstate_SBS_UE_int = reshape(sum(out,2),size(SBS_location,1),size(UE_location,1));
   Linkstate_SBS_UE = (Linkstate_SBS_UE_int==0);
   Linkstate_SBS_UE = Linkstate_SBS_UE';  % for dimension matching in later steps
    % For SBS-MBS
   MBS_loc_rep = repmat(MBS_location,size(SBS_location,1),1);
   SBS_loc_rep  = repelem(SBS_location,size(MBS_location,1),1);
   out= lineSegmentIntersect([SBS_loc_rep,MBS_loc_rep],[Block_endpoints_x(:,1),Block_endpoints_y(:,1),...
                  Block_endpoints_x(:,2),Block_endpoints_y(:,2)]);
   Linkstate_MBS_SBS_int = reshape(sum(out,2),size(MBS_location,1),size(SBS_location,1));
   Linkstate_MBS_SBS = (Linkstate_MBS_SBS_int==0);
   Linkstate_MBS_SBS = Linkstate_MBS_SBS';  % for dimension matching in later steps
   
   %% Step 3: Computing Association
   % Step 3-a: MBS
   X=UE_location';
   Y=MBS_location';
   a = sum(X.^2, 1);
   b = sum(Y.^2, 1);
   N1=size(X,2);
   N2=size(Y,2);
   P =  abs(a' * ones(1, N2)  + ones(N1, 1) * b - 2 * (X' * Y));
   P = sqrt(P); %to be consistent with pdist2
   SigPow_UE_MBS = P_m* L(P,Linkstate_MBS_UE);
   [Max_Power_UE_MBS, MBS_index] = max(SigPow_UE_MBS,[],2);
   % MBS_index holds the index of the candidate serving MBS
   % Step 3-b: SBS
   X=UE_location';
   Y=SBS_location';
   a = sum(X.^2, 1);
   b = sum(Y.^2, 1);
   N1=size(X,2);
   N2=size(Y,2);
   P =  abs(a' * ones(1, N2)  + ones(N1, 1) * b - 2 * (X' * Y));
   P = sqrt(P); %to be consistent with pdist2
   if is_shadowing
     Shadow_coeff_UE_MBS =   10.^(.1*(0 + Shadow_coefficient*randn(size(P))));
     SigPow_UE_SBS = P_s*Shadow_coeff_UE_MBS .* L(P);
   else
     SigPow_UE_SBS = P_s*L(P,Linkstate_SBS_UE);
   end
   [Max_Power_UE_SBS, SBS_index] = max(SigPow_UE_SBS,[],2);
   % SBS_index holds the index of the candidate serving SBS
   % Compute the association 
   Association_indicator = Max_Power_UE_MBS.*T_m>Max_Power_UE_SBS.*T_s;
   if Association_indicator(1) 
       internal_counter_m = internal_counter_m+1;
   end
   fprintf('\n Association Probability to MBS (A_m) = %f',internal_counter_m/count_sim);
end
fprintf('\n Association Probability to MBS (A_m) = %f',internal_counter_m/count_sim);