% Paramètres du modèle
ds = 0.1;   % Coefficient de diffusion pour s
di = 0.1;   % Coefficient de diffusion pour i
dr = 0.1;   % Coefficient de diffusion pour r
alpha = 0.3; % Paramètre alpha
beta2 = 0.1; % Paramètre beta2
beta3 = 0.1; % Paramètre beta3
gamma = 0.4; % Paramètre gamma
K = 100;      % Paramètre K
ksi = 0.5;  % Paramètre ksi

% Model parameters
N = 6*10^7; % Total population N = S + I + R
I0 = 100; % initial number of infected
T = 300; % period of 300 days
dt = 1/4; % time interval of 6 hours (1/4 of a day)
% fprintf('Value of parameter R0 is %.2f',N*beta/gamma)

% Calculate the model
[S,I,R] = sir_model(alpha, beta2, beta3, ksi,K, gamma,N,I0,T,dt);
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',2); grid on;
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R');

% Map
plot(I(1:(T/dt)-1),I(2:T/dt),"LineWidth",1,"Color",'r');
hold on; grid on;
plot(I(2),I(1),'ob','MarkerSize',4);
xlabel('Infected at time t'); ylabel('Infected at time t+1');
hold off;

function [S,I,R] = sir_model(alpha, beta2, beta3, ksi,K, gamma,N,I0,T,dt)
    % if delta = 0 we assume a model without immunity loss
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    for tt = 1:(T/dt)-1
        % Equations of the model
        dS = (-S(tt)*(K-S(tt))*(S(tt)-alpha)*I(tt) - gamma*S(tt)/(1-ksi*S(tt)) - beta2*S(tt)*R(tt)) * dt;
        dI = (S(tt)*(K-S(tt))*(S(tt)-alpha)*I(tt) - gamma*I(tt)/(1-ksi*I(tt)) - beta3*I(tt)*R(tt)) * dt;
        dR = (gamma*I(tt)/(1-ksi*I(tt)) +  gamma*S(tt)/(1-ksi*S(tt)) + beta2*S(tt)*R(tt) + beta3*I(tt)*R(tt)) * dt;
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
    end
end


