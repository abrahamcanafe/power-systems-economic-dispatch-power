%% EEE259_Final_Project_YALMIP.m
%
% Author: Abraham Canafe
%
% Date: April 14, 2022
%
%
% Switched from CPLEX to YALMIP.
%
%

clear, clc, close all
fullData = readmatrix('C:\Users\abrah\Dropbox (Personal)\Sac State\Spring 2022\EEE 259\Course Project.xlsx');
% cfs = cubic feet per second

%% Units Capacity Limits (MW)				
%           Pmin	Pmax	Start up Cost ($)	Duty Factors (cfs/MWh)
% Superior	0       20      1000                5
% High      0       75      2500                10
% Middle	0       45      2000                15
% Low       0       30      1500                25

%% Reservoir Storage Limits (ACF or Acre-feet per hour) *or cubic feet? (1 acre-foot = 43559.9 cubic-ft)
%               Min     Max     Initial	 End
% Roseville     8000	10000	9000	9000
% Sacramento	2500	3500	3000	3000
% Folsom        250     350     300     300
% El Dorado     2000	3000	2500	2500
% Rancho        400     600     500     500

%% LMP ($/MWh)
%           1    2    3     4    5    6     7     8     9      10    11      12      13      14      15      16      17      18      19      20      21      22      23      24
% Superior	$50  $45  $40 	$35  $40  $50 	$70   $90 	$50    $30 	 $10 	 $(10)	 $(20)	 $(20)	 $(5)	 $35 	 $40 	 $80 	 $100 	 $150 	 $125 	 $100 	 $80 	 $65 
% High      $52  $50  $50 	$30  $20  $80 	$90   $100 	$50    $30 	 $10 	 $5 	 $7 	 $10 	 $50 	 $55 	 $70 	 $100 	 $200 	 $250 	 $25 	 $200 	 $100 	 $75 
% Middle	$35  $30  $25 	$25  $30  $(10)	$(10) $(80)	$(70)  $(30) $(90)	 $(100)	 $(20)	 $(35)	 $(60)	 $(25)	 $10 	 $25 	 $50 	 $70 	 $72 	 $75 	 $60 	 $30 
% Low       $100 $125 $175 	$189 $50  $25 	$90   $130 	$49    $30 	 $25 	 $20 	 $25 	 $75 	 $85 	 $100 	 $140 	 $180 	 $250 	 $300 	 $325 	 $200 	 $150 	 $100 

%%
unitCapLimits = fullData(5:8,2:5);
reservoirStorage = fullData(13:17,2:5);
LMP = fullData(5:8,9:32);

Pmin = unitCapLimits(:,1); %+ 1e-10; % TESTING FOR BINARY VARIABLES
Pmax = unitCapLimits(:,2);
Rmin = Pmin;
Rmax = Pmax;
startupCosts = unitCapLimits(:,3);
rawDutyFactors = unitCapLimits(:,4);
minElevations = reservoirStorage(:,1);
maxElevations = reservoirStorage(:,2);
firstLastElevations = reservoirStorage(:,3); % first and last elevations are the same
% also, ignore the last elevation for Roseville

% Conversion constant
ACFH_PER_CFS = 1/12.1; % https://www.kylesconverter.com/flow/acre--feet-per-hour-to-cubic-feet-per-second
%%%%%% Magnitudes of natural flows (in ACF/h) %%%%%%
FlowA = 10 * ACFH_PER_CFS;
FlowB = 25 * ACFH_PER_CFS;
FlowC = 5 * ACFH_PER_CFS;
FlowD = 5 * ACFH_PER_CFS;
FlowE = 4.7 * ACFH_PER_CFS;
FlowF = 10 * ACFH_PER_CFS;
FlowG = 8.3 * ACFH_PER_CFS;
FlowH = 6.7 * ACFH_PER_CFS;
FlowI = 20 * ACFH_PER_CFS;
FlowJ = 9 * ACFH_PER_CFS;
FlowK = 30 * ACFH_PER_CFS;
FlowL = 400 * ACFH_PER_CFS;
FlowM = 8 * ACFH_PER_CFS;
FlowN = 10 * ACFH_PER_CFS;

% Assume that these net natural flows are inflows into the reservoirs.
totNatFlowL1 = FlowA + FlowB + FlowC;
totNatFlowL2 = -FlowD + FlowE + FlowF - FlowG;
totNatFlowL3 = FlowG + FlowH + FlowI - FlowJ;
totNatFlowL4 = FlowJ + FlowK - FlowL;
totNatFlowL5 = FlowL - FlowM - FlowN;

vecNatFlow = [totNatFlowL1;
              totNatFlowL2;
              totNatFlowL3;
              totNatFlowL4;
              totNatFlowL5];

% Variable counts (for CPLEX)
nUnits = 4; % flows in ACF/h
nPowers = nUnits; % power in MW
nReservoirs = 5;
nBinary = nUnits;
nSD = nUnits;
nSU = nUnits;
nReserve = nUnits; % new

nHours = 24;
maxAward = 1000;
totVars = (nUnits + nPowers + nReservoirs + nBinary + nSD + nSU + nReserve) * nHours;

%f = zeros(totVars,1);
%%%%%% Arrangement of variables for CPLEX (696) %%%%%%
% X1(t)     X2(t)       X3(t)       X4(t)                       [1 - 96]
% P1(t)     P2(t)       P3(t)       P4(t)                       [97 - 192]
% L1(t)     L2(t)       L3(t)       L4(t)       L5(t)           [193 - 312]
% u1(t)     u2(t)       u3(t)       u4(t)                       [313 - 408]    
% SD1(t)    SD2(t)      SD3(t)      SD4(t)                      [409 - 504]
% SU1(t)    SU2(t)      SU3(t)      SU4(t)                      [505 - 600]
% R1(t)     R2(t)       R3(t)       R4(t)                       [601 - 696]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pattern: X1(1), X2(1), X3(1), X4(1) -> X1(2), X2(2), X3(2), X4(2) -> ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPIN_PRICE = 5;
unwrappedLMP = reshape(LMP,[],1);

columnLMP = LMP';
startupCostVector = repmat(startupCosts, 1, nHours)'; % $
rawDutyFactorsVector = repmat(rawDutyFactors, 1, nHours)'; % cfs/MWh
minElevationsVector = repmat(minElevations, 1, nHours)'; % ACF
maxElevationsVector = repmat(maxElevations, 1, nHours)'; % ACF
firstLastElevationsVector = repmat(maxElevations, 1, nHours)'; % ACF
flowConversionVector = ACFH_PER_CFS*rawDutyFactorsVector; % ACF-h/cfs * cfs/MWh = ACF-h/MWh

%% CPLEX Variables
% Flows
X1 = sdpvar(nHours,1);
X2 = sdpvar(nHours,1);
X3 = sdpvar(nHours,1);
X4 = sdpvar(nHours,1);

% Power
P1 = sdpvar(nHours,1);
P2 = sdpvar(nHours,1);
P3 = sdpvar(nHours,1);
P4 = sdpvar(nHours,1);

% Reserve
R1 = sdpvar(nHours,1);
R2 = sdpvar(nHours,1);
R3 = sdpvar(nHours,1);
R4 = sdpvar(nHours,1);

% Reservoir Levels
L1 = sdpvar(nHours,1);
L2 = sdpvar(nHours,1);
L3 = sdpvar(nHours,1);
L4 = sdpvar(nHours,1);
L5 = sdpvar(nHours,1);

% Binary Variables
u1 = binvar(nHours,1);
u2 = binvar(nHours,1);
u3 = binvar(nHours,1);
u4 = binvar(nHours,1);

% Shutdown Binary Variables
SD1 = binvar(nHours,1);
SD2 = binvar(nHours,1);
SD3 = binvar(nHours,1);
SD4 = binvar(nHours,1);

% Startup Binary Variables
SU1 = binvar(nHours,1);
SU2 = binvar(nHours,1);
SU3 = binvar(nHours,1);
SU4 = binvar(nHours,1);

% Objective Function
Objective = -sum(columnLMP(:,1).*P1) - sum(columnLMP(:,2).*P2) - sum(columnLMP(:,3).*P3) - sum(columnLMP(:,4).*P4)...
    + sum(SPIN_PRICE*P1) + sum(SPIN_PRICE*P2) + sum(SPIN_PRICE*P3) + sum(SPIN_PRICE*P4)...
    + sum(SPIN_PRICE*R1) + sum(SPIN_PRICE*R2) + sum(SPIN_PRICE*R3) + sum(SPIN_PRICE*R4)...
    + sum(startupCostVector(:,1).*SU1) + sum(startupCostVector(:,2).*SU2)...
    + sum(startupCostVector(:,3).*SU3) + sum(startupCostVector(:,4).*SU4);

%% Flow = Power + Reserve
flow_equivalence = [flowConversionVector(1,1).*(P1(1) + R1(1)) == X1(1)];
for i = 2:nHours
    flow_equivalence = [flow_equivalence, flowConversionVector(i,1).*(P1(i) + R1(i)) == X1(i)];
end

flow_equivalence = [flow_equivalence, flowConversionVector(1,2).*(P2(1) + R2(1)) == X2(1)];
for i = 2:nHours
    flow_equivalence = [flow_equivalence, flowConversionVector(i,2).*(P2(i) + R2(i)) == X2(i)];
end

flow_equivalence = [flow_equivalence, flowConversionVector(1,3).*(P3(1) + R3(1)) == X3(1)];
for i = 2:nHours
    flow_equivalence = [flow_equivalence, flowConversionVector(i,3).*(P3(i) + R3(i)) == X3(i)];
end

flow_equivalence = [flow_equivalence, flowConversionVector(1,4).*(P4(1) + R4(1)) == X4(1)];
for i = 2:nHours
    flow_equivalence = [flow_equivalence, flowConversionVector(i,4).*(P4(i) + R4(i)) == X4(i)];
end

%% Power and Reserve Bounds
power_reserve_bounds = [Pmin(1,1)*u1(1) <= P1(1) <= Pmax(1,1)*u1(1)];
power_reserve_bounds = [power_reserve_bounds, Pmin(1,1)*u1(1) <= P1(1) + R1(1) <= Pmax(1,1)*u1(1)];
power_reserve_bounds = [power_reserve_bounds, 0 <= R1(1) <= Rmax(1,1)*u1(1)];
for i = 2:nHours
    power_reserve_bounds = [power_reserve_bounds, Pmin(1,1)*u1(i) <= P1(i) <= Pmax(1,1)*u1(i)];
    power_reserve_bounds = [power_reserve_bounds, Pmin(1,1)*u1(i) <= P1(i) + R1(i) <= Pmax(1,1)*u1(i)];
    power_reserve_bounds = [power_reserve_bounds, 0 <= R1(i) <= Rmax(1,1)*u1(i)];
end

power_reserve_bounds = [power_reserve_bounds, Pmin(2,1)*u2(1) <= P2(1) <= Pmax(2,1)*u2(1)];
power_reserve_bounds = [power_reserve_bounds, Pmin(2,1)*u2(1) <= P2(1) + R2(1) <= Pmax(2,1)*u2(1)];
power_reserve_bounds = [power_reserve_bounds, 0 <= R2(1) <= Rmax(2,1)*u2(1)];
for i = 2:nHours
    power_reserve_bounds = [power_reserve_bounds, Pmin(2,1)*u2(i) <= P2(i) <= Pmax(2,1)*u2(i)];
    power_reserve_bounds = [power_reserve_bounds, Pmin(2,1)*u2(i) <= P2(i) + R2(i) <= Pmax(2,1)*u2(i)];
    power_reserve_bounds = [power_reserve_bounds, 0 <= R2(i) <= Rmax(2,1)*u2(i)];
end

power_reserve_bounds = [power_reserve_bounds, Pmin(3,1)*u3(1) <= P3(1) <= Pmax(3,1)*u3(1)];
power_reserve_bounds = [power_reserve_bounds, Pmin(3,1)*u3(1) <= P3(1) + R3(1) <= Pmax(3,1)*u3(1)];
power_reserve_bounds = [power_reserve_bounds, 0 <= R3(1) <= Rmax(3,1)*u3(1)];
for i = 2:nHours
    power_reserve_bounds = [power_reserve_bounds, Pmin(3,1)*u3(i) <= P3(i) <= Pmax(3,1)*u3(i)];
    power_reserve_bounds = [power_reserve_bounds, Pmin(3,1)*u3(i) <= P3(i) + R3(i) <= Pmax(3,1)*u3(i)];
    power_reserve_bounds = [power_reserve_bounds, 0 <= R3(i) <= Rmax(3,1)*u3(i)];
end

power_reserve_bounds = [power_reserve_bounds, Pmin(4,1)*u4(1) <= P4(1) <= Pmax(4,1)*u4(1)];
power_reserve_bounds = [power_reserve_bounds, Pmin(4,1)*u4(1) <= P4(1) + R4(1) <= Pmax(4,1)*u4(1)];
power_reserve_bounds = [power_reserve_bounds, 0 <= R4(1) <= Rmax(4,1)*u4(1)];
for i = 2:nHours
    power_reserve_bounds = [power_reserve_bounds, Pmin(4,1)*u4(i) <= P4(i) <= Pmax(4,1)*u4(i)];
    power_reserve_bounds = [power_reserve_bounds, Pmin(4,1)*u4(i) <= P4(i) + R4(i) <= Pmax(4,1)*u4(i)];
    power_reserve_bounds = [power_reserve_bounds, 0 <= R4(i) <= Rmax(4,1)*u4(i)];
end

%% Awards limits
awards_limits = [sum(P1) + sum(P2) + sum(P3) + sum(P4) <= 1000];

%% Reservoir Flow Equations
reservoir_flow_equations = [L1(1) == firstLastElevations(1,1) - X1(1) + totNatFlowL1];
for i = 2:nHours
    reservoir_flow_equations = [reservoir_flow_equations, L1(i) == L1(i-1) - X1(i) + totNatFlowL1];
end

reservoir_flow_equations = [reservoir_flow_equations, L2(1) == firstLastElevations(2,1) + X1(1) - X2(1) + totNatFlowL2];
for i = 2:nHours
    reservoir_flow_equations = [reservoir_flow_equations, L2(i) == L2(i-1) + X1(i) - X2(i) + totNatFlowL2];
end

reservoir_flow_equations = [reservoir_flow_equations, L3(1) == firstLastElevations(3,1) + X2(1) - X3(1) + totNatFlowL3];
for i = 2:nHours
    reservoir_flow_equations = [reservoir_flow_equations, L3(i) == L3(i-1) + X2(i) - X3(i) + totNatFlowL3];
end

reservoir_flow_equations = [reservoir_flow_equations, L4(1) == firstLastElevations(4,1) + X3(1) + totNatFlowL4];
for i = 2:nHours
    reservoir_flow_equations = [reservoir_flow_equations, L4(i) == L4(i-1) + X3(i) + totNatFlowL4];
end

reservoir_flow_equations =  [reservoir_flow_equations, L5(1) == firstLastElevations(5,1) - X4(1) + totNatFlowL5];
for i = 2:nHours
    reservoir_flow_equations = [reservoir_flow_equations, L5(i) == L5(i-1) - X4(i) + totNatFlowL5];
end

%% Startup variable formulations
SU_SD_formulations = [u1(1) - SU1(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, -u1(1) - SD1(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, u1(1) - SU1(1) + SD1(1) <= 0];
for i = 2:nHours
    SU_SD_formulations = [SU_SD_formulations, u1(i) - u1(i-1) - SU1(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, -u1(i) + u1(i-1) - SD1(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, u1(i) - u1(i-1) - SU1(i) + SD1(i) <= 0];
end

SU_SD_formulations = [SU_SD_formulations, u2(1) - SU2(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, -u2(1) - SD2(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, u2(1) - SU2(1) + SD2(1) <= 0];
for i = 2:nHours
    SU_SD_formulations = [SU_SD_formulations, u2(i) - u2(i-1) - SU2(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, -u2(i) + u2(i-1) - SD2(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, u2(i) - u2(i-1) - SU2(i) + SD2(i) <= 0];
end

SU_SD_formulations = [SU_SD_formulations, u3(1) - SU3(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, -u3(1) - SD3(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, u3(1) - SU3(1) + SD3(1) <= 0];
for i = 2:nHours
    SU_SD_formulations = [SU_SD_formulations, u3(i) - u3(i-1) - SU3(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, -u3(i) + u3(i-1) - SD3(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, u3(i) - u3(i-1) - SU3(i) + SD3(i) <= 0];
end

SU_SD_formulations = [SU_SD_formulations, u4(1) - SU4(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, -u4(1) - SD4(1) <= 0];
SU_SD_formulations = [SU_SD_formulations, u4(1) - SU4(1) + SD4(1) <= 0];
for i = 2:nHours
    SU_SD_formulations = [SU_SD_formulations, u4(i) - u4(i-1) - SU4(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, -u4(i) + u4(i-1) - SD4(i) <= 0];
    SU_SD_formulations = [SU_SD_formulations, u4(i) - u4(i-1) - SU4(i) + SD4(i) <= 0];
end


%% Final reservoir values
%final_reservoir_values = [L1(nHours) == firstLastElevations(1,1)];

% Still not sure how to get all four reservoirs to have the desired storage
% levels at this point without getting an infeasible solution.
final_reservoir_values = [L2(nHours) == firstLastElevations(2,1)];
final_reservoir_values = [final_reservoir_values, L3(nHours) == firstLastElevations(3,1)];
%final_reservoir_values = [final_reservoir_values, L4(nHours) == firstLastElevations(4,1)];
final_reservoir_values = [final_reservoir_values, L5(nHours) == firstLastElevations(5,1)];

%%
options = sdpsettings('verbose',1,'solver','BNB');
sol = optimize([flow_equivalence, power_reserve_bounds, awards_limits,...
    reservoir_flow_equations, SU_SD_formulations, final_reservoir_values], Objective, options);
totalCost = value(Objective);
fprintf('\nNet cost: $%6.2f\n', totalCost);
fprintf('Profit: $%6.2f\n', -totalCost);

%result = solvesdp([flow_equivalence, power_reserve_bounds, awards_limits, reservoir_flow_equations],Objective,options)

results_X = value([X1 X2 X3 X4]);
results_P = value([P1 P2 P3 P4]);
results_R = value([R1 R2 R3 R4]);
results_u = value([u1 u2 u3 u4]);
results_L = value([L1 L2 L3 L4 L5]);
results_SU = value([SU1 SU2 SU3 SU4]);
results_SD = value([SD1 SD2 SD3 SD4]);

% Display final levels
fprintf('Final Roseville Reservoir Storage: %4.1f\n', results_L(nHours,1));
fprintf('Final Sacramento Reservoir Storage: %4.1f (Success)\n', results_L(nHours,2));
fprintf('Final Folsom Reservoir Storage: %4.1f (Success)\n', results_L(nHours,3));
fprintf('Final El Dorado Reservoir Storage: %4.1f (Fail)\n', results_L(nHours,4));
fprintf('Final Rancho Cordova Reservoir Storage: %4.1f (Success)\n', results_L(nHours,5));
fprintf('\nFlow, Power, Reserve, Binary, Reservoir, SU, and SD variables\n');
fprintf('are stored in ''results_X'', ''results_P'', ''results_X'', ''results_u'',\n');
fprintf('''results_L'', ''results_SU'', and ''results_SD''.\n');


