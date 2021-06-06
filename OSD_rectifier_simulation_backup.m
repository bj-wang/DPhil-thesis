% Simulation of the output voltage-frequency characteristics of OSD-based rectifiers
% Version 1.0, by Bingjun Wang, 16/Jul/2020, Department of Physics, University of Oxford, UK
% Version 2.0, by Bingjun Wang, 04/Apr/2021, Department of Physics, University of Oxford, UK
% Version 2.1, by Bingjun Wang, 04/May/2021, Department of Physics, University of Oxford, UK
% Version 2.2, by Bingjun Wang, 06/Jun/2021, Department of Physics, University of Oxford, UK
% Email: bingjun.wang@physics.ox.ac.uk or bingjun.wang1995@outlook.com
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License: https://creativecommons.org/licenses/by-nc-sa/4.0/.


%% Section 1: Initialisation

clearvars;              % clear all variables
tic


%% Section 2: Parameter Input

% Circuit parameter
V_m = 5;                % maximum value of the voltage, V

% Diode parameters
mu_raw = 10^(-4);       % mobility, cm^2 V^(-1) s^(-1). Default = 1E-4
d_raw = 100;            % layer thickness, nm. Default = 100
V_T = 0.5;              % turn-on voltage, V. Default = 0.5
J_r_raw = 10^(-4);      % leakage current density, A cm^(-2). Default = 1E-4
N_t = 5*10^23;          % trap density, m^(-3). Default = 5E23
m = 4;                  % trap depth parameter, m = dJ / dV. Default = 4


%% Section 3: Constants

ep_0 = 8.8542*10^(-12);     % vacuum permittivity, F m^(-1)
ep_r = 3.5;                 % relative permittivity
q = 1.6*10^(-19);           % elementry charge, C
N = 5*10^26;                % density of transport states, m^(-3)


%% Section 4: Simulation Settings

V_interval = 0.01;                    % interval of V_o, V
f_threshold = 1;                      % the minimum frequency at which the simulation finishes
output_name = 'standard';
output_path = 'file path';


%% Section 5: Conversion to SI Units and Deriving Other Parameters

mu = mu_raw * 10^(-4);
d = d_raw * 10^(-9);
J_r = J_r_raw * 10^4;
V_TC = q*d*d/(ep_0*ep_r) * (9/8 * N_t^m / N * ((m+1)/m)^m * ((m+1)/(2*m+1))^(m+1))^(1/(m-1));   % trap-controlled voltage, V
if V_TC < V_T
    V_TC = V_T;                                                                                 % make sure V_TC >= V_T
end
TCSCLC_pre = q*mu*N* ((ep_0*ep_r)/(q*N_t))^m *(m/(m+1))^m *((2*m+1)/(m+1))^(m+1);               % the prefactor in the TC-SCLC formula


%% Section 6: Simulation

M = zeros(V_m/V_interval*2,2);        % first column: V_o, V; second column: frequency, Hz. There are way more numbers of rows than needed to support the simulation in the next section.

for i=1:1:V_m/V_interval-1
    M(i,1) = V_m - i*V_interval;
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


%% Section 7: Achieving Desired Low-Frequency Response

% Find the smallest positive frequency obtained in the previous simulation
i = i + 1;
j = 1;
while M(j,2) < 0
    j = j + 1;
end
M(i,1) = M(j,1);
M(i,2) = M(j,2);

V_neg = M(j-1,1);       % the V_o of the frequency that is negative but closest to 0

while M(i,2) > f_threshold
    i = i + 1;
    M(i,1) = 0.5 * (M(i-1,1) + V_neg);      % use the bisection method to define V_o for the frequency calculation
    
    syms omega t Q_r J_1 J_2 Q_1 Q_2;
    Q_load = M(i,1) * ep_0 * ep_r / d;                                                      % charge per unit area consumed by the load per period
    Q_r = J_r*(asin(M(i,1)/V_m)/omega + 2*pi/omega - (pi-asin(M(i,1)/V_m))/omega);          % charge per unit area wasted by the diode's reverse leakage current per period
    
    if M(i,1) + V_TC >= V_m
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
    
    % If the obtained frequency is negative, set the corresponding V_o to be V_neg and repeat the calculation
    if M(i,2) < 0
        V_neg = M(i,1);
        i = i - 1;
    end
end


%% Section 8: Clean the Data and Plot the Figure

M(all(M == 0,2),:) = [];        % clear all zero rows
M((M(:,2) < 0),:) = [];         % clear all negative-frequency-rows
M(:,[1 2])=M(:,[2 1]);          % swap the dependent and the independent variables for the ease of plotting
M = sortrows(M);                % sort the rows in ascending order of frequency

semilogx(M(:,1),M(:,2))         % plot V_o vs f figure
xlabel('Frequency (Hz)');
ylabel('Vo (V)');


%% Section 9: Save the Data (for further analysis in Origin) and the Plot

output_table = array2table(M,'VariableNames',{'Frequency (Hz)','Vo (V)'});
file_path = pwd;
cd(output_path);
writetable(output_table,strcat(output_name,'.csv'));
print(gcf,strcat(output_name,'.tif'),'-dtiff', '-r300');
cd(file_path);
toc


%% Section 10: Sound notice

% The silumation typically takes more than 10 min using an Intel Core I7-3632QM CPU; sometimes it can take half an hour. Thus there is a sound to notify you when the simulation ends.
load train
sound(y,Fs)
