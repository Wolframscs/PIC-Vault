% 电极平衡案例-主函数
clc
clear all
%% parameter input
tic
global SOL_pos_ref;
global Epos_ref;
global SOL_neg_ref;
global Eneg_ref;
% global E_IR_pos;
% global E_IR_neg;
global SOC_exp;
global OCV_exp;
%global SOC;


%% data input
filepath_pos='.\LFP_OCV_5999_200.csv';
filepath_neg='.\QE_1.csv';
filepath_full_cell='.\FC_OCV_24599.xlsx';
% [SOL_pos_ref,Epos_ref]=textread(filepath_pos,'%f %f');
% [SOL_neg_ref,Eneg_ref]=textread(filepath_neg,'%f %f');
pos_data=csvread(filepath_pos);
[C,ia,ic]=unique(pos_data(:,1),'stable');
pos_data=pos_data(ia,:);
SOL_pos_ref=pos_data(:,1);
Epos_ref=pos_data(:,2);
neg_data=csvread(filepath_neg);
[C,ia,ic]=unique(neg_data(:,1),'stable');
neg_data=neg_data(ia,:);
SOL_neg_ref=neg_data(:,1);
Eneg_ref=neg_data(:,2);
% [data_full_SOC,data_full_V]=textread(filepath_full_cell,'%f %f');
% data_full=cat(2,data_full_SOC,data_full_V);

data_full=xlsread(filepath_full_cell);
[C,ia,ic]=unique(data_full(:,1),'stable');
data_updata=data_full(ia,:);
SOC_exp=data_updata(:,1);
OCV_exp=data_updata(:,2);
% SOC_exp=data_updata.data_updat(:,1);
% OCV_exp=data_updata.data_updat(:,2);
%% Problem Definition

CostFunction=@electrode_balance_example;        % Cost Function

nVar=6;            % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=[0.5 0.5 0.9 0.05 -0.02 -0.02];         % Lower Bound of Variables
VarMax= [0.7 0.75 1.02 0.08 0.02 0.02];         % Upper Bound of Variables


%% PSO Parameters

MaxIt=100;      % Maximum Number of Iterations

nPop=50;        % Population Size (Swarm Size)

% PSO Parameters
% w=1;            % Inertia Weight
% wdamp=0.99;     % Inertia Weight Damping Ratio
% c1=1.5;         % Personal Learning Coefficient
% c2=2;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w=w*wdamp;
    
end

BestSol = GlobalBest;

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
toc
plot_bestposition(BestSol.Position)


% 电极平衡目标函数
function [RMSE]=electrode_balance_example(x)
epss_pos=x(:,1);
epss_neg=x(:,2);
SOL_pos_0=x(:,3);
SOL_neg_0=x(:,4);
E_IR_pos=x(:,5);
E_IR_neg=x(:,6);
%%
Qcell=197.993; %Ah
S=479.5*112*69*2/1000000; %m^2
L_pos=7.52E-5; %m
L_neg=5.38E-5; %m
%epss_pos=0.58
%epss_neg=0.6
Fconst=96485; %C.mol-1
csmax_pos=22821; %mol.m-3
csmax_neg=31250; %mol.m-3
%SOL_pos_0=0.8
%SOL_neg_0=0.1
global SOL_pos_ref;
global Epos_ref;
global SOL_neg_ref;
global Eneg_ref;
% global E_IR_pos;
% global E_IR_neg;
global SOC_exp;
global OCV_exp;
% E_IR_pos=0; %ohm
% E_IR_neg=0; %ohm
%global SOC;
try
SOC=linspace(0,1,300);
SOL_pos=SOL_pos_0-SOC.*Qcell./((S*L_pos*epss_pos*csmax_pos*Fconst)/3600);
SOL_neg=SOL_neg_0+SOC.*Qcell./((S*L_neg*epss_neg*csmax_neg*Fconst)/3600);    
Ecell_sim=interp1(SOL_pos_ref,Epos_ref,SOL_pos,'pchip','extrap')+E_IR_pos...
    -(interp1(SOL_neg_ref,Eneg_ref,SOL_neg,'pchip','extrap')+E_IR_neg);
Ecell_exp=interp1(SOC_exp,OCV_exp,SOC,'pchip','extrap');
diff=(Ecell_exp-Ecell_sim).*1000;
RMSE=sqrt(mean(diff.*diff));
catch
    RMSE=666;
end
end


open('XP.fig');
obj=get(gca,'children');
x1=get(obj(1),'xdata')';
y1=get(obj(1),'ydata')';
y2=get(obj(2),'ydata')';


% 电极平衡绘图
function [] = plot_bestposition(x)
%PLOT_BESTPOSITION 此处显示有关此函数的摘要
%   此处显示详细说明
epss_pos=x(:,1);
epss_neg=x(:,2);
SOL_pos_0=x(:,3);
SOL_neg_0=x(:,4);
E_IR_pos=x(:,5);
E_IR_neg=x(:,6);
%%
Qcell=197.993; %Ah
S=479.5*112*69*2/1000000; %m^2
L_pos=7.52E-5; %m
L_neg=5.38E-5; %m
%epss_pos=0.58
%epss_neg=0.6
Fconst=96485; %C.mol-1
csmax_pos=22821; %mol.m-3
csmax_neg=31250; %mol.m-3
%SOL_pos_0=0.8
%SOL_neg_0=0.1
global SOL_pos_ref;
global Epos_ref;
global SOL_neg_ref;
global Eneg_ref;
global SOC_exp;
global OCV_exp;
%%
SOC=linspace(0,1,300);
SOL_pos=SOL_pos_0-SOC.*Qcell./((S*L_pos*epss_pos*csmax_pos*Fconst)/3600);
SOL_neg=SOL_neg_0+SOC.*Qcell./((S*L_neg*epss_neg*csmax_neg*Fconst)/3600);    
Ecell_sim=interp1(SOL_pos_ref,Epos_ref,SOL_pos,'pchip','extrap')+E_IR_pos...
    -(interp1(SOL_neg_ref,Eneg_ref,SOL_neg,'pchip','extrap')+E_IR_neg);
Ecell_exp=interp1(SOC_exp,OCV_exp,SOC,'pchip','extrap');
%%
figure;
hold on;grid on;
plot(SOC,Ecell_sim,'-b','LineWidth',2);
plot(SOC,Ecell_exp,'--or','LineWidth',1,'MarkerIndices',1:length(SOC)/50:length(SOC));
xlabel('SOC');ylabel('Voltage(V)');
legend('sim-data','exp-data','Location','SouthEast');
hold off;

end



