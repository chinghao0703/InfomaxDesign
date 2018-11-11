function dx=kinase_cascade_lin10(t,x,alpha, beta, C, R, Rt)

%% Variables
% alpha(i) is the rate constant for kinase reaction i
% beta(i)  is the rate constant for phosphatse reaction i
% C(i) is the total concentration of kinase i
% R is the receptor activation function which is evaluated at the time grid Rt = linspace(0,5,tspan);
% The function is taken to be R = exp(-lambda* Rt)


%% Here we need to interpolate the receptor activation function R
R = interp1(Rt,R,t); % Interpolate the data set (ft,f) at time t


% Equations are based on simple mass action kinetics 
% The exceptions are transcriptional induction of HSPs and Reporter, which use Hill function


dx(1) = alpha(1)*R.*(1-x(1)/C(1)) - beta(1)*x(1);

dx(2) = alpha(2)*x(1)*(1-x(2)/C(2))- beta(1)*x(2);

dx(3) = alpha(3)*x(2)*(1-x(3)/C(3))- beta(3)*x(3);

dx(4) = alpha(4)*x(3)*(1-x(4)/C(4))- beta(4)*x(4);

dx(5) = alpha(5)*x(4)*(1-x(5)/C(5))- beta(5)*x(5);

dx(6) = alpha(6)*x(5)*(1-x(6)/C(6))- beta(6)*x(6);

dx(7) = alpha(7)*x(6)*(1-x(7)/C(7))- beta(7)*x(7);

dx(8) = alpha(8)*x(7)*(1-x(8)/C(8))- beta(8)*x(8);

dx(9) = alpha(9)*x(8)*(1-x(9)/C(9))- beta(9)*x(9);

dx(10) = alpha(10)*x(9)*(1-x(10)/C(10))- beta(10)*x(10);



% matrix equation to be evaluated by ODE solver 
dx = [dx(1);dx(2);dx(3);dx(4);dx(5);dx(6);dx(7);dx(8);dx(9);dx(10)];

end