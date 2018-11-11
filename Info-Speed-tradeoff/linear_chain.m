clear all;
close all;
clc;

N = 10; % number of proteins in the linear chain
Nsamples = 100;
mu = linspace(-5, 0.1,100);
sigma = 0.1;
epsilon = 10^(-6);
q = [0.1, 0.2, 0.5, 0.7, 1.0];
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

            Pn= M*P1 + epsilon;
            M = M + epsilon;
            MI(iq,im,is) = sum((log2(M./Pn).*M)*P1);

        end

    end
    MImf(iq) = fm(q(iq), epsilon);
end

MImf(1,length(MImf)) = 0;

MImean = mean(MI, 3);
hold on;
plot(mu, MImean(1,:),'LineWidth', 3);
plot(mu, MImean(2,:),'LineWidth', 3);
plot(mu, MImean(3,:),'LineWidth', 3);
plot(mu, MImean(4,:),'LineWidth', 3);
plot(mu, MImean(5,:),'LineWidth', 3);
axis([-5 0 0 1]);
ax = gca;
legend('q=0.1','q=0.2', 'q=0.5', 'q=0.7', 'q=1.0');
set(gca,'FontSize',26)
xlabel('Mean binding affinity, $\mu$', 'Interpreter', 'latex');
ylabel('Mutual information, $I_{N=8}$', 'Interpreter','latex');


print('MI_v_theta','-depsc');
print('MI_v_theta','-dpdf','-bestfit')

%pd = histfit(MI);


function x = f(y)
% return the conditional probability (i.e. transfer matrix)

x = exp(-y)/ (1 + exp(-y)) ;


end

function x = fm(p, delta)
% mean-field approximation in the tight-binding limit
 
x= p*(1-delta)*log2((1-delta)/p) + (1-p)*delta*log2(delta/p) + p*delta*log2(delta/(1-p)) + (1-p)*(1-delta)*log2((1-delta)/(1-p));

end


