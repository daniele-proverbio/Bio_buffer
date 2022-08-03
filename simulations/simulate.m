%% Optimise EWS for detection of Fold CT
% Looking for te best combination of summary statistis to detect critical
% transitions driven by fold bifurcations, in a biological toy model
% (autoactivating feedback loop motif for gene regulation)

%% Author
% Daniele Proverbio, 04/08/2021
% daniele.proverbio@outlook.com / @uni.lu
% University of Luxembourg


%% Prepare env
clear; close all; clc;


%% Initialize

simu = 5;

enne = [2,3,4,5,8];
p_critical = [1.788,1.737,1.62,1.524,1.344];  % obtained from studies of bifurcation diagram

N_Exp = 200;          %repeated experiments

time_start = 1;     % Start time of simulations 
N= 10000;            % Time points
dt = 0.01;           % Time step 
tic = 0;

val1 = 0.1;               %0.1 normal % Basal expression (constant, within accepted range), normally = 0.1 ; otherwise = 0  
val2 = p_critical(simu)+0.35:-0.005:p_critical(simu)+0.15;     % Max Production (control parameter)  
noise = 0.02;             % Noise level (diffusion term) -> lower than basin height

dist = val2 - p_critical(simu)+0.01;

sol = zeros((N*dt+1)/dt,length(val2),N_Exp);

%% SDE simulator
% Euler Maruyama scheme
for experiment =1 : N_Exp                        %repeated experiments
    
    x_in = 1.3; %+ (2.25-1.3)*rand;             % Initial condition (on upper branch) .  (initial conditions: x_in = [2.25, 1.3])
    
    for m = 1 : length(val2)   
        tic = 0;        
        for p = time_start-1+dt : dt : N*dt+1
            tic = tic + 1;          
            if tic == 1
                sol(tic,m,experiment) = x_in;
            else
                f = protein_production(p-dt,sol(tic-1,m,experiment),val1,val2(m),enne(simu));
                sol(tic,m,experiment) = sol(tic-1,m,experiment) + f * dt + noise*sqrt(dt)*randn;
            end            
        end
    end
end


save(['multiple_exps_ct_wn_n', num2str(enne(simu)) ,'.mat'],'sol', '-v7.3')

%% Equation
% Simulate simple equation for autocatalytic protein production (deterministic part) (Sharma,
% 2015)

function dxdt = protein_production(t,x,val1,val2,enne)

K = val1;    
c = val2;    

dxdt = K + (c*x^enne)/(1+x^enne) - x;

end
