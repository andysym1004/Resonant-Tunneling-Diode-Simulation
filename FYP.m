
%electrical constants
e = 1.60217662e-19;  % electron charge (C)
hbar = 1.0545718e-34; % reduced Planck's constant (Js)
m0 = 9.10938356e-31;   % electron mass (kg)
k_B = 1.38e-23; % Boltzmann constant (J/K)


%input parameters
Lw = 5e-9; % quantum well width (m)
Lb = 2e-9; % barrier width (m)
Le = 2e-8; % emitter width (m)
Lc = 2e-8; % collector width (m)
ND_GaAs= 1e18;  % doping concentration of GaAs emmiter and collector (cm^-3)
T = 300; % Temperature (K)
x = 0.3; % mole fraction composition of AlxGa(1-x)As


%material properties for GaAs (emmitter, collector, well)
m_GaAs = 0.067 * m0; % approximate effective mass of GaAs
Nc_GaAs = 4.7e17 * (T/300)^(3/2);  % density of states in conduction band in cm^-3
Eg0_GaAs = 1.523;  % Bandgap energy at 0K in eV
alpha_GaAs = 5.4e-4;
beta_GaAs = 204;
Eg_GaAs= Eg0_GaAs - (alpha_GaAs * T^2) / (T + beta_GaAs);
Ec_GaAS = Eg_GaAs;  % Assuming conduction band minimum at Eg
Ef_GaAs = Ec_GaAS + k_B/e * T * log(ND_GaAs / Nc_GaAs);  


%material properties for AlxGa(1-x)As (barrier)
m_AlGaAs = (0.067 + 0.083*x) * m0; % effective mass for GaAs (kg)
Nc_AlGaAs = 4.82e15 * (T/300)^(3/2);  % density of states in conduction band in cm^-3
Eg0_AlGaAs = 1.523 + 1.23*x;  
alpha_AlGaAs = (5.5 + 3.35*x + 50*x^2) * 10^-4;
beta_AlGaAs = 198 + 88*x + 4300*x^2;
Eg_AlGaAs = Eg0_AlGaAs - (alpha_AlGaAs * T^2) / (T + beta_AlGaAs);
Ec_AlGaAS = Eg_AlGaAs;  % Assuming conduction band minimum at Eg

%grid parameterisation:
N = 1000; % no. of discretised points
dn = 1;
V_bias= 0.0; % applied biased voltage
n = linspace(0, Le+ 2*Lb + Lw + Lc, N); % grid range
ne = sum(n >= 0 & n < Le);
nc = N - ne + 1;
nbw = sum(n >= Le & n <= Le+ 2*Lb + Lw);

%effective mass
m = zeros(size(n));
for i = 1:length(n)
    if n(i) < Le
        m(i) = m_GaAs;  % emitter
    elseif n(i) >= Le && n(i) < Le + Lb
        m(i) = m_AlGaAs; % left barrier 
    elseif n(i) >= Le + Lb && n(i) <= Le + Lb +Lw
        m(i) = m_GaAs; % well 
    elseif n(i) > Le + Lb +Lw && n(i) <= Le + 2*Lb +Lw
        m(i) = m_AlGaAs; %right barrier Vb +
    else
        m(i) = m_GaAs;  % collector
    end
end

%potential profile
V = zeros(size(N));
for i = 1:length(n)
    if n(i) < Le
        V(i) = Ec_GaAS + V_bias;  % emitter
    elseif n(i) >= Le && n(i) < Le + Lb
        V(i) = Ec_AlGaAS + V_bias * (N - (i-ne)*N/nbw)/N; % left barrier 
    elseif n(i) >= Le + Lb && n(i) <= Le + Lb +Lw
        V(i) = Ec_GaAS + V_bias * (N - (i-ne)*N/nbw)/N; % well 
    elseif n(i) > Le + Lb +Lw && n(i) <= Le + 2*Lb +Lw
        V(i) = Ec_AlGaAS + V_bias * (N - (i-ne)*N/nbw)/N;  %right barrier Vb +
    else
        V(i) = Ec_GaAS;  % collector
    end
end
% 
% %Hamiltonian operator
% H0 = zeros(N);
% H1 = zeros(N);
% H2 = zeros(N);
% 
% 
% for i = 1:N
%     H0(i, i) = (hbar^2 / (m(i) * dn^2)) + V(i);  % Main diagonal elements
%     if i < N
%         H1(i, i+1) = -hbar^2 / (2 * m(i) * dn^2);  % Upper diagonal elements
%         H2(i+1, i) = -hbar^2 / (2 * m(i+1) * dn^2);  % Lower diagonal elements
%     end
% end
% [fn,En] = eig(H);
% 
% Er = diag(En);



%J = current_density(V, Er, e, Gamma, m_star, k_B, T, E_F)


%plot energy band profile
Efplot = zeros(size(n)); % data for plotting Ef for emitter and collector region
for i = 1:length(n)
    if n(i) < Le
       Efplot(i) = Ef_GaAs + V_bias;
    elseif n(i) > Le + 2*Lb +Lw
        Efplot(i) = Ef_GaAs;
    else 
        Efplot(i) = 0;
    end
end

figure;
plot(n, V);
hold on;
 plot (n(1:ne),Efplot(1:ne),n(nc:N),Efplot(nc:N));

 hold off;
xlabel('Position (m)');
ylabel('Energy (eV)');
title('Energy Profile of RTD under applied bias');
grid on;






%solve s


%current density 
function J = current_density(V, Er, e, Gamma, m_star, k_B, T, E_F)
    T = @(E, V) (Gamma / 2)^2 ./ ((E - (Er - e * V / 2)).^2 + (Gamma / 2)^2);
    J_integrand = @(E) T(E, V) .* log((1 + exp((E_F - E) / (k_B * T))) ./ ...
        (1 + exp((E_F - E - e * V) / (k_B * T))));
    J_integral = (e * m_star * k_B * T) ./ (2 * pi^2 * h_bar^3) * integral(@(E) J_integrand(E), 0, inf);
end





%functions
% %transmission probability
% function T = transmission(E, V, Er, e, Gamma)
% 
%     % Calculate the transmission probability T(E, V)
%     T = ((Gamma / 2)^2) / (E - (Er - e .* V / 2))^2 + (Gamma ./ 2)^2;
% 
% end
% 
% 


