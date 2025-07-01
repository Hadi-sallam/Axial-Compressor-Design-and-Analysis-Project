%% cod e initialization
clearvars
clc
%% inputs
pressure_ratio = 3;    % At SLS conditions.
m_dot = 0.5;            % kg/s
N = 30000;             % RPM
U_over_rm = 2*pi*N/60;
h_over_c = 2;
M_rel_max = 0.75; % which is also Mrel_1 as the compressor slows down the flow
DF_max=0.5;
zeta = 0.75;
T1 = 288;              % sls conditions
P1 = 1.01325*10^5;
gamma = 1.4;
R_air = 287;
Cp = 1004.5;

phi_range  = linspace(0.4,0.6,10); %setting the values up for iterations later
psi_range  = linspace(0.2,0.4,10);
R_range    = linspace(0.5,0.8,10);
DF_range   = [0.4,0.5];
h_c_range  = [2,3.5];


%% optimizing for efficiency
eff=0;
% the loop keeps iterating to calculate the efficiency and stores the maximum effeciency calculated
for phi=phi_range                         % and the values that achieve it.
    for psi=psi_range                     % efficiency(phi,psi,R) is simply a snippet of this code up to the effeciency calculation.
        for R=R_range                     % to make the code run a bit faster.
            f=efficiency(phi,psi,R);
            if double(f) > double(eff)
               eff=double(f);
               best_comb=[phi,psi,R];
            end           
        end
    end
end
%% Velocity triangles
phi=best_comb(1);                          % the phi, psi and R values calculated above are used to run the calculations needed.
psi=best_comb(2);
R=best_comb(3);
beta1  = atand(((R+0.5*psi)/phi));
beta2  = atand((R-0.5.*psi)/phi);
alpha1 = atand((1-R-0.5.*psi)/phi);
alpha2 = atand((1-R+0.5.*psi)/phi);
alpha3 = alpha1;  

Cx = M_rel_max * cosd(beta1) * sqrt(gamma * R_air * T1);
U = Cx*(tand(alpha1)+tand(beta1));


C1 = Cx/cosd(alpha1);
C2 = Cx/cosd(alpha2);
C_theta1 = Cx*tand(alpha1);
C_theta2 = Cx*tand(alpha2);

W1 = Cx/cosd(beta1);
W2 = Cx/cosd(beta2);
W_theta1 = U-C_theta1;
W_theta2 = U-C_theta2;

%% temperature calculation
Tt1 = T1+((C1^2)/(2*Cp));
Ttrel1=T1+(W1^2/(2*Cp));

Ttrel2=Ttrel1;
T2=Ttrel2-(W2^2/(2*Cp));
Tt2 = T2+((C2^2)/(2*Cp));

Tt3 = Tt2;             %The flow across the stator is adiabatic
T3=Tt3-((C1^2)/(2*Cp));
Ttrel3=T3+(W1^2/(2*Cp));

temp_ratio = Tt3/Tt1;

%% pressure calculation
Pt1=P1*((Tt1/T1)^(gamma/(gamma-1)));
Ptrel1=Pt1*((Ttrel1/Tt1)^(gamma/(gamma-1)));

R_solidity=(W_theta1-W_theta2)/((DF_max-1+(W2/W1))*2*W1)
S_solidity=(C_theta2-C_theta1)/((DF_max-1+(C1/C2))*2*C2)

LC_r=0.032;   % From figure of DF
LC_s=0.027;   % From figure of DF
w_r=(LC_r*2*R_solidity)/cosd(beta2);
w_s=(LC_s*2*S_solidity)/cosd(alpha1);

Ptrel2=Ptrel1-(w_r*(Ptrel1-P1));
Pt2=Ptrel2*(Tt2/Ttrel2)^(gamma/(gamma-1));
P2=Pt1*(T2/Tt1)^(gamma/(gamma-1));
Pt3=Pt2-(w_s*(Pt2-P2));
P3=Pt3*(T3/Tt3)^(gamma/(gamma-1));
Ptrel3 = P3*(Ttrel3/T3)^(gamma/(gamma-1));
Pr_stage1=Pt3/Pt1;
%% work and effieciency calculation
w_sp = U*(C_theta2-C_theta1);
Power_stage = m_dot*w_sp;
eta=(((Pr_stage1^((gamma-1)/gamma))-1))/(temp_ratio-1); % for a single stage
%% compressor design
delta_T_Stage=Tt3-Tt1;
N_stage=(Tt1*0.85*((pressure_ratio^((gamma-1)/gamma))-1))/delta_T_Stage;
N_stage=ceil(N_stage)

%Areas
A1 = m_dot/((P1/(R_air*T1))*Cx);
A2 = m_dot/((P2/(R_air*T2))*Cx);
A3 = m_dot/((P3/(R_air*T3))*Cx);


rm=U/U_over_rm;
An = 4*pi*rm^2*(1-zeta)/(1+zeta)
h = An/(2*pi*rm);
Chord = h/h_over_c;
N_Blades = (2*pi*rm*R_solidity)/Chord;
N_Blades = ceil(N_Blades)
compressor_Length = (cosd(beta2)+cosd(alpha2))*Chord*N_stage

%% Temperatures for different stages

Tt_vec=zeros(N_stage*2+1,1);
Tt_vec(1)=Tt1;
Tt_vec(2)=Tt2;
Tt_vec(3)=Tt3;

T_vec=zeros(N_stage*2+1,1);
T_vec(1)=T1;
T_vec(2)=T2;
T_vec(3)=T3;

Ttr_vec=zeros(N_stage*2+1,1);
Ttr_vec(1)=Ttrel1;
Ttr_vec(2)=Ttrel2;
Ttr_vec(3)=Ttrel3;

for j=4:size(Tt_vec)
    if mod(j,2) ~= 0            % when j is even we are calculating across the rotor and when it's odd
       Tt_vec(j)=Tt_vec(j-1);   % we are calculating across the stator.
       T_vec(j) = Tt_vec(j) - C1^2/(2*Cp);
       Ttr_vec(j) = T_vec(j) + W1^2/(2*Cp);
    else
       Tt_vec(j)=Tt_vec(j-1) + delta_T_Stage;
       T_vec(j) = Tt_vec(j) - C2^2/(2*Cp);
       Ttr_vec(j) = T_vec(j) + W2^2/(2*Cp);
    end
end
Tt_vec
Ttr_vec
T_vec
%% pressures

Pt_vec=zeros(N_stage*2+1,1);
Pt_vec(1)=Pt1;
Pt_vec(2)=Pt2;
Pt_vec(3)=Pt3;

P_vec=zeros(N_stage*2+1,1);
P_vec(1)=P1;
P_vec(2)=P2;
P_vec(3)=P3;

Ptr_vec=zeros(N_stage*2+1,1);
Ptr_vec(1)=Ptrel1;
Ptr_vec(2)=Ptrel2;
Ptr_vec(3)=Ptrel3;

for j=4:size(Pt_vec)
    if mod(j,2) ~=0             % same as the temperature vecotrs above
    Pt_vec(j)=Pt_vec(j-1)-(w_s*(Pt_vec(j-1)-P_vec(j-1)));
    P_vec(j)=Pt_vec(j)*(T_vec(j)/Tt_vec(j))^(gamma/(gamma-1));
    Ptr_vec(j) = Pt_vec(j)*(Ttr_vec(j)/Tt_vec(j))^(gamma/(gamma-1));
    else
    Ptr_vec(j)=Ptr_vec(j-1)-(w_r*(Ptr_vec(j-1)-P_vec(j-1)));
    P_vec(j)=Ptr_vec(j)*(T_vec(j)/Ttr_vec(j))^(gamma/(gamma-1));
    Pt_vec(j) = P_vec(j)*(Tt_vec(j)/T_vec(j))^(gamma/(gamma-1));
    end
end
Pt_vec
Ptr_vec
P_vec
%% total power ana total efficiency
w_sp_total = Cp*(Tt_vec(13) - Tt_vec(1))
Total_Power = m_dot*w_sp_total
tau_total = Tt_vec(13) / Tt_vec(1);
pi_total = Pt_vec(13) / Pt_vec(1);
efficiency_total =(((pi_total^((gamma-1)/gamma))-1))/(tau_total-1)


%% Mach numbers

M1 = C1/sqrt(gamma*R_air*T1);
M = zeros(N_stage*2+1,1);
Mrel = zeros(N_stage*2+1,1);
M(1) = M1;
Mrel(1) = M_rel_max;

for j = 2:size(M)
    if mod(j,2) ~= 0
        M(j) = C1/sqrt(gamma*R_air*T_vec(j));
        Mrel(j) =  W1/sqrt(gamma*R_air*T_vec(j));
    else
        M(j) = C2/sqrt(gamma*R_air*T_vec(j));
        Mrel(j) =W2/sqrt(gamma*R_air*T_vec(j));
    end
end
M
Mrel
%% blades angle

% rotor
n=-0.28;
i_0=1;
a_c=0.5;
sigma_r=0.3369;
epsilon_r=beta1-beta2;
m_r=0.23*(2*a_c)^2+beta2/500;
syms beta1_dash beta2_dash
EQN1 = 1+n*(beta1_dash-beta2_dash) + beta1_dash == beta1;
EQN2 = m_r/sqrt(sigma_r)*(beta1_dash-beta2_dash) + beta2_dash == beta2;
sol_r = solve([EQN1,EQN2],[beta1_dash,beta2_dash]);
beta1_dash=double(sol_r.beta1_dash)
beta2_dash=double(sol_r.beta2_dash)

I_Rotor = beta1 - beta1_dash;                       % Rotor Incidence angle
Dev_angle_Rotor = beta2- beta2_dash ;               % Rotor Deviation angle
Theta_Rotor = beta1_dash - beta2_dash;              % Rotor Camber angle 
Stagger_Rotor = (beta1_dash+beta2_dash)/2;          % Rotor Stagger angle
AOA_Rotor = beta1-Stagger_Rotor;                    % Angle of Attack

% stator
sigma_s=0.3396;
espsilon_s= alpha2-alpha1;
m_s=0.23*(2*a_c)^2 + alpha1/500;
syms alpha2_dash alpha3_dash
EQN3= 1 + n*(alpha2_dash-alpha3_dash) + alpha2_dash == alpha2;
EQN4= m_s/sqrt(sigma_s)*(alpha2_dash-alpha3_dash) +alpha3_dash == alpha1;
sol_s = solve([EQN3,EQN4],[alpha2_dash,alpha3_dash]);
alpha2_dash=double(sol_s.alpha2_dash)
alpha3_dash=double(sol_s.alpha3_dash)

I_Stator = alpha2-alpha2_dash     ;             % Stator Incidence Angle
Dev_angle_Stator = alpha3-alpha3_dash  ;        % Stator Deviation Angle
Theta_Stator = alpha2_dash-alpha3_dash ;        % Stator Camber Angle 
Stagger_Stator=(alpha2_dash+alpha3_dash)/2;     % Stator stagger angle
AOA_Stator=alpha2-Stagger_Stator;                    % Stator angle of attack