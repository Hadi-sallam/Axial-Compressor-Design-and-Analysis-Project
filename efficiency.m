function eta=efficiency(phi,psi,R)
pressure_ratio = 3;    % At SLS conditions.
m_dot = 0.5;            % kg/s
N = 30000;             % RPM
U_over_rm = 2*pi*N/60;
M_rel_max = 0.75;
T1 = 288; 
P1 = 1.01325*10^5;
gamma = 1.4;
R_air = 287;
Cp = 1004.5;
DF_max=0.5;

beta1  = atand(((R+0.5*psi)/phi));
beta2  = atand((R-0.5.*psi)/phi);
alpha1 = atand((1-R-0.5.*psi)/phi);
alpha2 = atand((1-R+0.5.*psi)/phi);

% syms x
% lhs = x / cos(beta1) / sqrt(gamma * R_air * (Tt1 - x^2 / (2 * Cp * (cos(alpha1))^2)));
% rhs = M_rel_max;
% 
% 
% eqn_squared = lhs^2 == rhs^2;
% 
% x=solve(eqn_squared,x,"Real",true);
% cx=double(x);
% Cx=max(cx);
W1= M_rel_max*sqrt(gamma*R_air*T1);
Cx = W1*cosd(beta1);
W2 = Cx/cosd(beta2);
C1 = Cx/cosd(alpha1);
C2 = Cx/cosd(alpha2);
U = Cx*(tand(alpha1)+tand(beta1));
C_theta1 = Cx*tand(alpha1);
C_theta2 = Cx*tand(alpha2);
W_theta1 = U-C_theta1;
W_theta2 = U-C_theta2;

Tt1 = T1+((C1^2)/(2*Cp));
Ttrel1=T1+(W1^2/(2*Cp));
Ttrel2=Ttrel1;
T2=Ttrel2-(W2^2/(2*Cp));
Tt2 = T2+((C2^2)/(2*Cp));
Tt3 = Tt2;             %The flow across the stator is adiabatic
T3=Tt3-((C1^2)/(2*Cp));
temp_ratio = Tt3/Tt1;

Pt1=P1*((Tt1/T1)^(gamma/(gamma-1)));
Ptrel1=Pt1*((Ttrel1/Tt1)^(gamma/(gamma-1)));

R_solidity=(W_theta1-W_theta2)/((DF_max-1+(W2/W1))*2*W1);
S_solidity=(C_theta2-C_theta1)/((DF_max-1+(C1/C2))*2*C2);

LC_r=0.032;   % From figure of DF
LC_s=0.027;   % From figure of DF
w_r=(LC_r*2*R_solidity)/cosd(beta2);
w_s=(LC_s*2*S_solidity)/cosd(alpha1);

Ptrel2=Ptrel1-(w_r*(Ptrel1-P1));
Pt2=Ptrel2*(Tt2/Ttrel2)^(gamma/(gamma-1));
P2=Pt1*(T2/Tt1)^(gamma/(gamma-1));
Pt3=Pt2-(w_s*(Pt2-P2));
P3=Pt3*(T3/Tt3)^(gamma/(gamma-1));
Pr_stage1=Pt3/Pt1;

w_sp = U*(C_theta2-C_theta1);
W_dot = m_dot*w_sp;

eta=double((((Pr_stage1^((gamma-1)/gamma))-1))/(temp_ratio-1));

% delta_Tt=Tt3-Tt1;
% N_stage=(-Tt1+(((pressure_ratio^((gamma-1)/gamma)-1)/eta_max)+1)*Tt1)/delta_Tt; % We take N=4
% 
% A1 = m_dot/((P1/(R_air*T1))*Cx);
% A2 = m_dot/((P2/(R_air*T2))*Cx);
% A3 = m_dot/((P3/(R_air*T3))*Cx);
% area_ratio=A2/A1 ;
% 
% rm=U/U_over_rm;
% h1=A1/(2*pi*rm);
% c=h1/h_c_range(1);
% r_hub=(2*rm-h1)/2;
% r_tip=h1+r_hub;
% Z=ceil(2*pi*rm*R_solidity)/c;
end