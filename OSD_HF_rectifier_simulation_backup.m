% Bingjun Wang, Department of Physics, University of Oxford, UK, April 2021
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License: https://creativecommons.org/licenses/by-nc-sa/3.0/.

%% Initialisation
clearvars;              % clear all variables
tic

%% Circuit parameters
V_m = 5;                % maximum value of the voltage, V

%% Diode parameters
mu_raw = 10^(-4);       % mobility, cm^2 V^(-1) s^(-1). Default = 1E-4
d_raw = 100;            % layer thickness, nm. Default = 100
V_T = 0.5;              % turn-on voltage, V. Default = 0.5
J_r_raw = 10^(-4);      % leakage current density, A cm^(-2). Default = 1E-4
N_t = 5*10^23;          % trap density, m^(-3). Default = 5E23
m = 4;                  % trap depth parameter, m = dJ / dV. Default = 4

%% Constants
ep_0 = 8.8542*10^(-12);     % vacuum permittivity, F m^(-1)
ep_r = 3.5;                 % relative permittivity
q = 1.6*10^(-19);           % elementry charge, C
N = 5*10^26;                % density of transport states, m^(-3)

%% Conversion to SI units and deriving other parameters
mu = mu_raw * 10^(-4);
d = d_raw * 10^(-9);
J_r = J_r_raw * 10^4;
V_TC = q*d*d/(ep_0*ep_r) * (9/8 * N_t^m / N * ((m+1)/m)^m * ((m+1)/(2*m+1))^(m+1))^(1/(m-1));   % trap-controlled voltage, V
if V_TC < V_T
    V_TC = V_T;                                                                                 % make sure V_TC >= V_T
end
TCSCLC_pre = q*mu*N* ((ep_0*ep_r)/(q*N_t))^m *(m/(m+1))^m *((2*m+1)/(m+1))^(m+1);               % the prefactor in the TC-SCLC formula

%% Simulation
interval = 0.02;                    % interval of V_o, V
M = zeros(V_m/interval-1,2);        % first column: V_o, V; second column: frequency, Hz

for i=1:1:V_m/interval-1
    M(i,1) = V_m-i*interval;
    syms omega t Q_r J_1 J_2 Q_1 Q_2;
    Q_load = M(i,1) * ep_0 * ep_r / d;                                                      % charge per unit area consumed by the load per period
    Q_r = J_r*(asin(M(i,1)/V_m)/omega + 2*pi/omega - (pi-asin(M(i,1)/V_m))/omega);          % charge per unit area wasted by the diode's reverse leakage current per period
    
    if M(i,1) + V_T >= V_m
        M(i,2) = -1;                                                                        % no corresponding frequency
                
    elseif M(i,1) + V_TC >= V_m
        J_1 = TCSCLC_pre *(V_m*sin(omega*t)-M(i,1)).^(m+1) / (d^(2*m+1));                   % trap-controlled SCLC
        Q_1 = 2 * int(J_1,t,asin((M(i,1)+V_T)/V_m)/omega,pi/(2*omega));                     % charge per unit area provided by the diode in the TC-SCLC regime per period
        M(i,2) = double(solve(Q_1 - Q_load - Q_r,omega))/(2*pi);
        
    else
        J_1 = TCSCLC_pre *(V_m*sin(omega*t)-M(i,1)).^(m+1) / (d^(2*m+1));                   % trap-controlled SCLC
        J_2 = 9/8 * ep_0*ep_r*mu * (V_m*sin(omega*t)-M(i,1)).^2 / (d^3);                    % trap-free SCLC
        Q_1 = 2 * int(J_1,t,asin((M(i,1)+V_T)/V_m)/omega,asin((M(i,1)+V_TC)/V_m)/omega);    % charge per unit area provided by the diode in the TC-SCLC regime per period
        Q_2 = 2 * int(J_2,t,asin((M(i,1)+V_TC)/V_m)/omega,pi/(2*omega));                    % charge per unit area provided by the diode in the TF-SCLC regime per period
        M(i,2) = double(solve(Q_1 + Q_2 - Q_load - Q_r,omega))/(2*pi);
    end
end

semilogx(M(:,2),M(:,1))         % plot V_o vs f figure
toc
