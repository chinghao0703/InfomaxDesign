clear all;
close all;
load('tradeoff_all_data.mat');




figure(1);
subplot(121)
loc = find(tau_mean(10,:)<0,1)-2;
yyaxis left
for i = [10]
    plot(theta(1:loc),I05(1:loc) , 'LineWidth', 3);
    hold on;
end
set(gca,'FontSize',26);
ylabel('Mutual information, $I$ (bits)', 'Interpreter', 'latex');
xlabel('Mean binding affinity, $\mu$ ', 'Interpreter', 'latex');

yyaxis right
for i = [10]
    plot(theta(1:loc), 1000./tau_mean(i,1:loc) ,'LineWidth', 3);
    hold on;
end
ylabel('speed, $\tau^{-1}$ (ms)', 'Interpreter', 'latex');


subplot(122)
plot(1000./tau_mean(i,1:loc),I05(1:loc) , 'k--','LineWidth', 3);
set(gca,'FontSize',26);
col = get(gca,'colororder');
ylabel('Mutual information, $I$ (bits)', 'Interpreter', 'latex', 'Color', col(1,:));
xlabel('Speed, $\tau^{-1}$ (ms)', 'Interpreter', 'latex', 'Color', col(2,:));

set(gcf, 'Position', [ 270         417        1174         480])

print('tradeoff_all','-depsc');
print('tradeoff_all','-djpeg');