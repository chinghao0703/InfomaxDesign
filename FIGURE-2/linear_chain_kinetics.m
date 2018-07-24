clear all;
clc;


load('MI_v_mu_lin_q05.mat');

delta = 10^(-3);
Tspan = 400;
Nsamples = 100;
C =0.2*(1+ zeros(1,10)); % total concentration of kinase
alpha = 1+ zeros(1,10); % phosphorylation rate 
%theta = linspace(-3,-0.5, 50);
theta = [-3];
sigma = 0.1;
activefraction = 0.1;

x0 = activefraction*C; % initial concentration
%x0 = [activefraction*C(1) 0 0 0];
tspan = [0 Tspan];


Rt = linspace(0,Tspan,2500);
t_half = round(0.5*length(Rt));
R = linspace(1, 1, length(Rt)); % 0.1 phase-shift to ensure at the off state
R(1,t_half:length(Rt)) = 0;
%R = exp(-lambda*Rt);

tau = zeros(10,length(theta),Nsamples);
tau2 = zeros(10,length(theta),Nsamples);
xmax = zeros(length(theta),Nsamples);



%% Note the ode function syntax 
% function dx=kinase_cascade(t,x,alpha, beta, C, R, Rt)


for im = 1:length(theta)
    
    mu = theta(im);
    
    for is = 1:Nsamples
        
        thetavec = mu + sigma*randn(1, 10);
        beta = exp(thetavec)./C.*alpha; % dephosphorylation rate
        
        opts = odeset('RelTol',1e-4,'AbsTol',1e-7);
        [t,x] = ode45(@(t,x)kinase_cascade_lin10(t,x,alpha, beta, C, R, Rt), tspan, x0, opts);
        
        
        for ix = 1:10
            [xm ,ti] = max(x(:,ix));
            ti2 = round(length(t)/2);
            t_temp = find(abs(x(:,ix))<= delta, 1);
            
            if  isempty(t_temp)
                t_temp = length(t);
            end
            tau(ix,im,is)= t(t_temp) - t(ti2);
            
        end
        
        xmax(im, is) = max(max(x));
        
        
    end
    
end


tau_mean = mean(tau, 3);
%% Below plot figure


figure(1);
%epsilon_b = fliplr(-theta); 
loc = find(tau_mean(10,:)<0,1)-2;
yyaxis left
for i = [10]
    plot(1000./tau_mean(i,1:loc),I05(1:loc) , 'LineWidth', 2);
    hold on;
end
set(gca,'FontSize',20);
ylabel('mutual information, $I$ (bits)', 'Interpreter', 'latex');
xlabel('speed, $\tau^{-1}$ (ms)', 'Interpreter', 'latex');

yyaxis right
for i = [10]
    plot(1000./tau_mean(i,1:loc), theta(1:loc) ,'LineWidth', 2);
    hold on;
end
ylabel('mean binding affinity, $\mu$ ', 'Interpreter', 'latex');
%legend('boxoff')
%set(h,...
%    'Position',[0.351785714285714 0.608928571428572 0.1125 0.308333333333333]);


%print('tradeoff_all','-depsc');
%print('tradeoff_all','-dpdf','-bestfit')
%print('tradeoff_all','-djpeg')

%save('tradeoff_all_data.mat');


%c = get(gca,'colororder');
figure(2)
for ind = [1, length(x0)]
    plot(t, (x(:,ind)/max(xmax)), 'LineWidth',2);
    hold on;
end
plot(Rt, R, 'k--' );
axis([0 tspan(2) 0 1.1]);
legend('i=1','i=10');
set(gca,'FontSize',18)
%title('Strongly activated CFFL')
xlabel('time, t','Interpreter', 'latex')
ylabel('accitivy of kinase $i$','Interpreter', 'latex')





