%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Slava Ermolaev
% Date: 06/16/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all %#ok<*CLALL>
close all
clc

%% Problem 2a

%  Initiate constants

%  Pressure
P = 101.325; % (kPa)

%  Error difference
error = 1e-12;

%  Liquid compositions vector
%  Variables: x_w - water, x_e - ethanol, x_a - acetone, x_b - 1-butanol
x_w = 0.21; x_e = 0.19; x_a = 0.27; x_b = 0.33;

%  Ideal conditions
% x_w = 0; x_e = 0.5; x_a = 0; x_b = 0.5;

% Test case compositions
% x_w = 0; x_e = 0; x_a = 0.47; x_b = 0.53;

X = [x_w; x_e; x_a; x_b];
N = length(X);

%  Constants for the Antoine Equation
A_vec = [16.3872; 16.8958; 14.3145; 15.3144];
B_vec = [3885.7; 3795.17; 2756.22; 3212.43];
C_vec = [230.17; 230.918; 228.06; 182.739];

%  Calculate vector of saturated temperatures by calling function CalcTsats
T_sats = CalcTsats(P, A_vec, B_vec, C_vec);

%  Calculate estimate of temperature
T = sum(X.*T_sats);

%  Initiate zero for initial temperature variable
T0 = 0;

%  Start Newton's method iterations
while abs(T-T0) > 1e-12

    %  Set new temperature as previous
    T0 = T;

    %  Calculate next temperature
    T = TbubbleNewton(X,T0,P,A_vec,B_vec,C_vec);

end

%  Calculate vapor bubble composition
Y = X.*CalcPsats(A_vec, B_vec, C_vec, T)/P;

%  Store solutions in SOL_q2a_T and SOL_q2a_y
SOL_q2a_T = T;
SOL_q2a_y = Y.';

%  Print bubble temperature
fprintf('2a:\nT_bubble = %.2f C\n', SOL_q2a_T)
fprintf('Composition: ')
SOL_q2a_y

%% Problem 2b

%  Set liquid composition
Y = [x_w; x_e; x_a; x_b];

%  Calculate estimate of temperature
T = sum(X.*T_sats);

%  Initiate zero for initial temperature variable
T0 = 0;

%  Start Newton's method iterations
while abs(T-T0) > 1e-12

    %  Set new temperature as previous
    T0 = T;

    %  Calculate next temperature
    T = TdewNewton(Y,T0,P,A_vec,B_vec,C_vec);

end

%  Calculate liquid dew composition
X = Y*P./CalcPsats(A_vec, B_vec, C_vec, T);

%  Store solutions in SOL_q2b_T and SOL_q2b_x
SOL_q2b_T = T;
SOL_q2b_x = X.';

%  Print dew temperature
fprintf('2b:\nT_dew = %.2f C\n', SOL_q2b_T)
fprintf('Composition: ')
SOL_q2b_x

%% Problem 3a - calculate bubble T

X = [x_w; x_e; x_a; x_b];

%  Calculate estimate of temperature
T = sum(X.*T_sats);

%  Calculate estimate of saturated pressures
P_sats = CalcPsats(A_vec,B_vec,C_vec,T);

%  Initialize all initial Phi to 1
Phi = ones(N,1);

%  Calculate gamma by UNIFAC model
gamma = UNIFAC(X,T);
gamma = gamma.';

%  Set acetone [P_sats(3)] as species j
Pj_sat = P/sum(X.*gamma./Phi.*P_sats/P_sats(3));

%  Calculate T using Pj_sat
T = B_vec(3)/(A_vec(3)-log(Pj_sat)) - C_vec(3);

%  Set initial temperature as zero
T0 = 0;

%  Iterate until temperature difference is less than error
while abs(T - T0) > error
    
    %  Set new temperature as initial
    T0 = T;

    %  Calculate new P_sats with new T
    P_sats = CalcPsats(A_vec,B_vec,C_vec,T0);

    %  Calculate vapor composition
    Y = (X.*gamma.*P_sats)./(Phi*P);

    %  Calculate new Phi with new T
    Phi = VirialCalc(Y,P,T0,P_sats);
    Phi = Phi.';

    %  Calculate new gamma with new T
    gamma = UNIFAC(X,T0);
    gamma = gamma.';

    %  Calculate new Pj_sat
    Pj_sat = P/sum(X.*gamma./Phi.*P_sats/P_sats(3));

    %  Calculate new T
    T = B_vec(3)/(A_vec(3)-log(Pj_sat)) - C_vec(3);

end

%  Calculate bubble point vapor composition
Y = X.*gamma.*P_sats./(Phi*P);

%  Save solutions to solution variables
SOL_q3_bublT = T;
SOL_q3_y = Y.';

%  Print bubble temperature
fprintf('3a:\nT_bubble = %.2f C\n', SOL_q3_bublT)
fprintf('Composition: ')
SOL_q3_y

%% Problem 3a - calculate dew T

%  Set vapor composition as liquid composition from last part
Y = X;

%  Set Phi and gamma to 1
Phi = ones(N,1);
gamma = ones(N,1);

%  Calculate T_sats
T_sats = CalcTsats(P,A_vec,B_vec,C_vec);

%  Calculate temperature estimate
T = sum(Y.*T_sats);

%  Calculate P_sats using estimate T
P_sats = CalcPsats(A_vec,B_vec,C_vec,T);

%  Calculate Pj_sat using acetone (j=3) as species j
Pj_sat = P*sum(Y.*Phi./gamma*P_sats(3)./P_sats);

%  Calculate T using Pj_sat
T = B_vec(3)/(A_vec(3)-log(Pj_sat)) - C_vec(3);

%  Set initial temperature to zero
T0 = 0;

while abs(T-T0) > error

    %  Set previous temprature as initial
    T0 = T;

    %  Calculate P_sats
    P_sats = CalcPsats(A_vec,B_vec,C_vec,T0);

    %  Calculate Phi
    Phi = VirialCalc(Y,P,T0,P_sats);
    Phi = Phi.';

    %  Set initial gamma to zero
    gamma0 = zeros(N,1);

    while abs(max(abs(gamma)) - max(abs(gamma0))) > error

        %  Set previous gamma as initial
        gamma0 = gamma;

        %  Calculate liquid composition
        X = Y.*Phi*P./gamma0./P_sats;

        %  Normalize X
        X = X/sum(X);

        %  Calculate new gamma
        gamma = UNIFAC(X,T0);
        gamma = gamma.';

    end

    %  Calculate Pj_sat using acetone (j=3) as species j
    Pj_sat = P*sum(Y.*Phi./gamma*P_sats(3)./P_sats);

    %  Calculate new temperature
    T = B_vec(3)/(A_vec(3)-log(Pj_sat)) - C_vec(3);

end
    
%  Store dew temperature in solution variable
SOL_q3_dewT = T;
%  Store liquid composition in solution variable
SOL_q3_x = X.';

%  Print results
fprintf('3a:\nT_dew = %.2f C\n', SOL_q3_dewT)
fprintf('Composition: ')
SOL_q3_x

%% Problem 4

%  Set temperature to 360 K (86.85 C)
T = 360 - 273.15;

%  Test case temperature
%T = 100;

%  Set mixture composition
Z = [x_w; x_e; x_a; x_b];

%  Calculating dew pressure

%  Set vapor composition as mixture composition
Y = Z;

%  Set Phi and gamma to one
Phi_dew = ones(N,1);
gamma_dew = ones(N,1);

%  Calculate P_sats
P_sats = CalcPsats(A_vec,B_vec,C_vec,T);

%  Calculate dew pressure
P_dew = 1/sum(Y.*Phi_dew./gamma_dew./P_sats);

%  Calculate liquid composition
X = Y.*Phi_dew*P_dew./gamma_dew./P_sats;

%  Calculate gamma
gamma_dew = UNIFAC(X,T);
gamma_dew = gamma_dew.';

%  Calculate dew pressure
P_dew = 1/sum(Y.*Phi_dew./gamma_dew./P_sats);

%  Set initial P_dew to zero
P_dew0 = 0;

while abs(P_dew - P_dew0) > error

    P_dew0 = P_dew;

    %  Calculate Phi
    Phi_dew = VirialCalc(Y,P_dew0,T,P_sats);
    Phi_dew = Phi_dew.';

    gamma_dew0 = 0;

    while sum(abs(gamma_dew - gamma_dew0) > error) ~= 0
        
        gamma_dew0 = gamma_dew;

        X = Y.*Phi_dew*P_dew0./gamma_dew0./P_sats;
        gamma_dew = UNIFAC(X,T);
        gamma_dew = gamma_dew.';
    end

    P_dew = 1/sum(Y.*Phi_dew./gamma_dew./P_sats);

end

%  Calculating bubble pressure

X = Z;
Phi_bubble = ones(N,1);
gamma_bubble = UNIFAC(X,T);
gamma_bubble = gamma_bubble.';
P_bubble = sum(X.*gamma_bubble.*P_sats./Phi_bubble);
P_bubble0 = 0;

while abs(P_bubble - P_bubble0) > error
    
    P_bubble0 = P_bubble;
    Y = X.*gamma_bubble.*P_sats./Phi_bubble/P_bubble0;
    Phi_bubble = VirialCalc(Y,P_bubble0,T,P_sats);
    Phi_bubble = Phi_bubble.';

    P_bubble = sum(X.*gamma_bubble.*P_sats./Phi_bubble);

end

if P_dew < P && P < P_bubble

    X = Z;
    Y = Z;

    gamma = (P-P_dew)*(gamma_bubble-gamma_dew)/(P_bubble-P_dew)+gamma_dew;
    Phi = (P-P_dew)*(Phi_bubble-Phi_dew)/(P_bubble-P_dew)+Phi_dew;

    V_dew = 1;
    V_bubble = 0;
    V = (P-P_dew)*(V_bubble-V_dew)/(P_bubble-P_dew)+V_dew;

    V0 = 0;
    X0 = zeros(N,1);
    Y0 = zeros(N,1);

    while sum(abs(X-X0)>error)>0 || sum(abs(Y-Y0)>error)>0 || abs(V-V0)>error

        X0 = X;
        Y0 = Y;
        V0 = V;

        K = gamma.*P_sats./Phi/P;
        F = sum(Z.*(K-1)./(1+V0.*(K-1)));
        dFdV = -sum(Z.*(K-1).^2./(1+V0.*(K-1)).^2);

        V = V0 - F./dFdV;

        while abs(V-V0)>error

            V0 = V;
            F = sum(Z.*(K-1)./(1+V0.*(K-1)));
            dFdV = -sum(Z.*(K-1).^2./(1+V0.*(K-1)).^2);
            V = V0 - F./dFdV;
        end

        X = Z./(1+V.*(K-1));
        Y = K.*X;
        gamma = UNIFAC(X,T);
        gamma = gamma.';
        Phi = VirialCalc(Y,P,T,P_sats);
        Phi = Phi.';
    end

    SOL_q4_y = Y.';
    SOL_q4_x = X.';
    SOL_q4_V = V.';
    SOL_q4_L = (1 - V).';

    fprintf('4:\n')
    fprintf('Vapor phase composition:')
    SOL_q4_y
    fprintf('Liquid phase composition:')
    SOL_q4_x
    fprintf('Vapor phase mole fraction: %.4f\n', SOL_q4_V)
    fprintf('Liquid phase mole fraction composition: %.4f\n', SOL_q4_L)
    
else

    fprintf('Mixture is not within the two phase region. Cannot compute.\n')
end

%% Functions

function T_sats = CalcTsats(P, A_vec, B_vec, C_vec)

%  This function calculates vector of saturated temperatures using Antoine
%  Equation and provided vectors of A, B, and C constants and pressure

T_sats = B_vec./(A_vec-log(P))-C_vec;

end

function P_sats = CalcPsats(A_vec, B_vec, C_vec, T)

%  This function calculates P_sats using Antoine Equation

P_sats = exp(A_vec - B_vec./(T+C_vec));

end

function dP_sats = CalcdPsats(A_vec, B_vec, C_vec, T)

%  This function calculates vector of derivatives of saturated pressures
%  with respect to temperature

dP_sats = exp(A_vec-B_vec./(T+C_vec)).*(B_vec./(T+C_vec).^2);

end

function F2a = f2a(X, P, P_sats)

%  This function returns the sum(X*P_sats) - P for Newton's method purposes

F2a = sum(X.*P_sats) - P;

end

function dF2a = df2a(X, dP_sats)

%  This function returns the sum(X*dP_sats) for Newton's method purposes

dF2a = sum(X.*dP_sats);

end

function T = TbubbleNewton(X,T,P,A_vec,B_vec,C_vec)

%  This function solves one iteration of Newton's method to calculate
%  temperature of bubble point

P_sats = CalcPsats(A_vec, B_vec, C_vec, T);
dP_sats = CalcdPsats(A_vec,B_vec,C_vec,T);
F = f2a(X,P,P_sats);
dF = df2a(X,dP_sats);

T = T - F/dF;

end

function F2b = f2b(Y, P, P_sats)

%  This function returns the sum(Y./P_sats)-1/P for Newton's method purposes
F2b = sum(Y./P_sats) - 1/P;

end

function dF2a = df2b(Y, P_sats, dP_sats)

%  This function returns the sum(-Y./P_sats.^2)*dP_sats for Newton's method purposes

dF2a = sum(-Y./P_sats.^2.*dP_sats);

end

function T = TdewNewton(Y,T,P,A_vec,B_vec,C_vec)

%  This function solves one interation of Newton's method to calculate
%  temperature of dew point

P_sats = CalcPsats(A_vec, B_vec, C_vec, T);
dP_sats = CalcdPsats(A_vec,B_vec,C_vec,T);
F = f2b(Y,P,P_sats);
dF= df2b(Y,P_sats,dP_sats);

T = T - F/dF;

end

function gamma = UNIFAC(X,T)

%  This function calculates activity coefficient gamma using UNIFAC model
%  It accepts composition vector of liquid phase X and temperature in Celcius

%  Transpose X to make it 1xN
X = X.';

%  Number of species
N = length(X);

%  Species indices (i)
%  1 - Water (H2O), 2 - Ethanol, 3 - Acetone, 4 - 1-butanol

%  Subgroups indices (k)
%  1 - CH3, 2 - CH2, 3 - OH, 4 - H2O, 5 - CH3CO

%  Number of subgroups
numSbgrps = 5;

%  Convert Celcius to Kelvin
T = T + 273.15;

%  Number of subgroups in each species, v(k,i)
%  Rows (k) - subgroups. Columns (i) - species
v = [0 1 1 1;
     0 1 0 3;
     0 1 0 1;
     1 0 0 0;
     0 0 1 0];

%  R(k) parameters
R = [0.9011 0.6744 1.0000 0.92 1.6724];

%  Q(k) parameters
Q = [0.848 0.540 1.200 1.4 1.448];

%  Energy interaction parameters matrix a(i,j)
a = [0.0000 0.0000  986.50 1318.0  476.40;
     0.0000 0.0000  986.50 1318.0  476.40;
     156.40 156.40  0.0000 353.50  84.000;
     300.00 300.00 -229.10 0.0000 -195.40; 
     26.76  26.76   164.50 472.50  0.0000];

%  Calculate r(i) for each species
r = R*v;

%  Calculate q(i) for each species
q = Q*v;

%  Preallocate e(numSbgrps,N)
e = zeros(numSbgrps,N);

%  Calculate e(k,i)
for i = 1:N
    for k = 1:numSbgrps
        e(k,i) = v(k,i)*Q(k)/q(i);
    end
end

%  Allocate theta(numSbgrps)
theta = zeros(1,numSbgrps);

%  Calculate theta(numSbgrps)
for k = 1:numSbgrps
    theta(k) = sum(X.*q.*e(k,:))./sum(X.*q);
end

%  Preallocate tau
tau = zeros(numSbgrps,numSbgrps);

%  Calculate tau(numSubgrps,numSbgrps);
for m = 1:numSbgrps
    tau(m,:) = exp(-a(m,:)/T);
end

%  Preallocate beta
beta = zeros(N,numSbgrps);

%  Calculate beta(N,numSbgrps)
for i = 1:N
    for k = 1:numSbgrps
        beta(i,k) = sum(e(:,i).*tau(:,k));
    end
end

%  Calculate s(numSbgrps)
s = theta*tau;

%  Calculate J(N)
J = r/sum(r.*X);

%  Calculate L(N)
L = q/sum(q.*X);

%  Calculate gammaC
gammaC = exp(1-J+log(J)-5*q.*(1-J./L+log(J./L)));

%  Preallocate gammaR
gammaR = zeros(1,N);

%  Calculate gammaR
for i = 1:N
    gammaR_sum = 0;

    for k = 1:numSbgrps
        gammaR_sum = gammaR_sum + (theta(k)*beta(i,k)/s(k)-e(k,i)*log(beta(i,k)/s(k)));
    end

    gammaR(i) = exp(q(i)*(1-gammaR_sum));

end

%  Calculate gamma
gamma = gammaC.*gammaR;

end

function Phi = VirialCalc(Y,P,T,P_sats)

%  This function calculates fugacity coefficient, Phi, using virial
%  correlation for species:
%  1 - water, 2 - ethanol, 3 - acetone, 4 - 1-butanol

%  Number of species
N = length(Y);

%  Initialize parameters
T = T + 273.15; % Conver C to K
R = 8.314;
omegas = [0.345 0.645 0.307 0.594];
Tc = [647.1 513.9 508.2 563.1];
Zc = [0.229 0.240 0.233 0.260];
Vc = [55.9 167 209 275];

%  Preallocate matrices for parameters
Tc_ij = zeros(N,N);
omegas_ij = zeros(N,N);
Zc_ij = zeros(N,N);
Vc_ij = zeros(N,N);
delta_ij = zeros(N,N);
Phi_i = zeros(1,N);

%  Calculate parameters
for i = 1:N
    for j = 1:N

        Tc_ij(i,j) = sqrt(Tc(i)*Tc(j));
        omegas_ij(i,j) = (omegas(i)+omegas(j))/2;
        Zc_ij(i,j) = (Zc(i)+Zc(j))/2;
        Vc_ij(i,j) = 0.001*((Vc(i)^(1/3)+Vc(j)^(1/3))/2)^3;

    end
end

Tr_ij = T./Tc_ij;
Pc_ij = Zc_ij*R.*Tc_ij./Vc_ij;

%  Calculate virial coefficients
B0_ij = -0.422./(Tr_ij).^(1.6) + 0.083;
B1_ij = -0.172./(Tr_ij).^(4.2) + 0.139;
Bhat_ij = B0_ij + omegas_ij.*B1_ij;
B_ij = R*Tc_ij.*Bhat_ij./Pc_ij;

%  Calculate deltas
for i = 1:N
    for j = 1:N
        if i == j % If diagonal element
            delta_ij(i,j) = 0;
        else
            delta_ij(i,j) = 2*B_ij(i,j) - B_ij(i,i) - B_ij(j,j);
        end
    end
end

%  Calculate Phis
for i = 1:N
    for j = 1:N
        for k = 1:N
            Phi_i(i) = Phi_i(i) + Y(j)*Y(k)*(2*delta_ij(j,i)-delta_ij(j,k));
        end
    end

    Phi_i(i) = exp((B_ij(i,i)*(P-P_sats(i))+0.5*P*Phi_i(i))/(R*T));
end

Phi = Phi_i;

end