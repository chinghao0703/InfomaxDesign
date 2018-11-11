clear all;
close all;
clc;

N = 10; % number of proteins in the linear chain
Nsamples = 1000;
mu = linspace(-5, -0.1,50);
sigma = 0.1;
epsilon = 10^(-6);
q = linspace(0, 1, 51);
MImf = zeros(1,5);




MI = zeros(length(q),length(mu),Nsamples);


for iq = 1:length(q)
    
    P1 = [q(iq); 1-q(iq)];
    
    for im = 1:length(mu)

        for is = 1: Nsamples

            Ebind = mu(im) + sigma*randn(1, N-1);
            M = [1, 0; 0, 1]; % P(xn|x1)

            for i = 1: N-1
                Mtemp = [f(Ebind(i)), 0; 1-f(Ebind(i)), 1];
                M = Mtemp*M;
            end

            Pn= M*P1;
            M(1,2) = epsilon ;
            MI(iq,im,is) = sum((log2(M./Pn).*M)*P1);

        end

    end
    MImf(iq) = fm(q(iq), epsilon);
end

MImf(1,length(MImf)) = 0;

MImean = mean(MI, 3);
[val, loc] = max(MImean, [], 1);
qloc = q(loc);


load('opt_q.mat');

figure(1);
[Q, MU]= meshgrid(q,mu);
%revMU= flipud(MU);
%murev = fliplr(-mu);
hold on;
pcolor(Q, MU, MImean');
scatter(qloc, mu, 'w', 'filled', 'LineWidth',2);
scatter(qstarmean, mu, 'r');
shading flat;
colorbar;
colormap parula;
axis([0.05 1 -5.0 -0.1]);
xlabel('input probability, $q$', 'Interpreter', 'latex');
ylabel('mean binding affinity, $\mu$', 'Interpreter', 'latex');
h = gca;
h_colorbar.FontSize = 18;
set(gca,'Fontsize', 20);
title('mutual information, $I_{N=8}$', 'Interpreter','latex');






print('MIpd_lin_q_theta','-depsc');
print('MIpd_lin_q_theta','-djpeg');

%pd = histfit(MI);

function x = f(y)
% return the conditional probability (i.e. transfer matrix)

x = exp(-y)/ (1 + exp(-y)) ;


end

function x = fm(p, delta)
% mean-field approximation in the tight-binding limit
 
x= p*(1-delta)*log2((1-delta)/p) + (1-p)*delta*log2(delta/p) + p*delta*log2(delta/(1-p)) + (1-p)*(1-delta)*log2((1-delta)/(1-p));

end
