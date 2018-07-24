clear all;
close all;
clc;

theta = linspace(-5,-1,11); % binding affinity
eta = linspace(-5,5.5,10); % x-talk level

delta = 10^(-6);

q1 = linspace(0.01,0.99,20);
q2 = linspace(0.01,0.99,20);
q3 = linspace(0.01,0.99,20);

[Q1, Q2, Q3] = meshgrid(q1,q2,q3);

Q1r = reshape(permute(Q1, [2 1 3]),[],1);
Q2r = reshape(permute(Q2, [2 1 3]),[],1);
Q3r = reshape(permute(Q3, [2 1 3]),[],1);



row2Delete = find(1-Q1-Q2-Q3 < 0);
Q1r(row2Delete, :) = [];
Q2r(row2Delete, :) = [];
Q3r(row2Delete, :) = [];


Q1r = [Q1r ; 0.25];
Q2r = [Q2r ; 0.25];
Q3r = [Q3r ; 0.25];
Qall = [Q1r Q2r Q3r];


C = Q1.*(1-Q1-Q2-Q3)-Q2.*Q3;  % dim: (q2,q1,q3)
Cr = reshape(permute(C, [2 1 3]),[],1);
Cr(row2Delete, :) = [];
Cr = [Cr ; 0];
    
    
MI = zeros(length(Q1r),length(theta), length(eta));
MI_label = 'q_theta_eta';




for iq = 1:length(Q1r)
    
    p12 = [Q1r(iq); Q2r(iq); Q3r(iq); 1-Q1r(iq)-Q2r(iq)-Q3r(iq)];
    
    for it = 1:length(theta)
        
        for ie = 1:length(eta)
            
            eb1 = theta(it);
            eb2 = theta(it) + eta(ie);
            M = [1 (1-f(eb1))*(1-f(eb2)) (1-f(eb2))*(1-f(eb1)) (1-g(eb1,eb2))^2;
                0 f(eb1)*(1-f(eb2)) f(eb2)*(1-f(eb1)) g(eb1,eb2)*(1-g(eb1,eb2));
                0 (1-f(eb1))*f(eb2) (1-f(eb2))*f(eb1) (1-g(eb1, eb2))*g(eb1, eb2);
                0 f(eb1)*f(eb2) f(eb2)*f(eb1) g(eb1,eb2)^2];
            p34 = M*p12;
            Md = M;
            Md(:,1) = [1 delta delta delta]';
            MI(iq, it, ie) = sum((log2(abs(Md./p34)).*M)*p12);
            
            
            
            
            
        end
        
    end
end








lpos = find(Cr >0);
lneg = find(Cr <0);
[v,l]=min(abs(Cr));
[vcmin,lcmin]=min(Cr);
[vcmax,lcmax]=max(Cr);


[v_minTB, l_minTB] = min(MI(:,1,1)-MI(:,1,length(eta)));
[v_minWB, l_minWB] = min(MI(:,length(theta),1)-MI(:,length(theta),length(eta)));

[v_maxTB, l_maxTB] = max(MI(:,1,1)-MI(:,1,length(eta)));
[v_maxWB, l_maxWB] = max(MI(:,length(theta),1)-MI(:,length(theta),length(eta)));



%{
%% figure 1-- most negative c 

hfig = figure(1); 
plot(eta, reshape(MI(lmin,1,:), [1, length(eta)]),'LineWidth', 2);% tight-binding
hold on;
plot(eta, reshape(MI(lmin, length(theta),:), [1, length(eta)]), 'LineWidth', 2');% weak-binding
hold on;


y1 = 2+zeros(1,length(eta));
H1=area(eta,y1);
hold on
ide=eta>0;
H=area(eta(ide),y1(ide));
set(H(1),'FaceColor',[1 0.5 0]);
alpha(0.125)
set(gca,'FontSize',18)
legend({'tight-binding', 'weak-binding'}, 'Location', 'northwest');
xlabel('cross-talk level, $\eta$', 'Interpreter', 'latex');
ylabel('mutual information, $I$', 'Interpreter','latex');
title('$c<0$','Interpreter', 'latex', 'Fontsize', 15)
txt1 = 'no cross-talk';
txt2 = 'cross-talk';
t= text(3, 1.5, txt1, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
t= text(-8, 1.5, txt2, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
hold off;

print('MI_xtalk_most_negC' ,'-depsc');
print('MI_xtalk_most_negC' ,'-djpeg');



%% figure 2-- most positive c 

hfig2 = figure(2);  
plot(eta, reshape(MI(lmax,1,:), [1, length(eta)]),'LineWidth', 2);% tight-binding
hold on;
plot(eta, reshape(MI(lmax, length(theta),:), [1, length(eta)]), 'LineWidth', 2');% weak-binding
hold on;



y1 = 2+zeros(1,length(eta));
H1=area(eta,y1);
hold on
ide=eta>0;
H=area(eta(ide),y1(ide));
set(H(1),'FaceColor',[1 0.5 0]);
alpha(0.125)
set(gca,'FontSize',18)
legend({'tight-binding', 'weak-binding'}, 'Location', 'northwest');
xlabel('cross-talk level, $\eta$', 'Interpreter', 'latex');
ylabel('mutual information, $I$', 'Interpreter','latex');
title('$c>0$','Interpreter', 'latex', 'Fontsize', 15)
txt1 = 'no cross-talk';
txt2 = 'cross-talk';
t= text(3, 0.3, txt1, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
t= text(-8, 0.3, txt2, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
hold off;

print('MI_xtalk_most_posC' ,'-depsc');
print('MI_xtalk_most_posC' ,'-djpeg');



%% figure 3 -- c = 0  

hfig3 = figure(3); 
[v,l]=min(abs(Cr));
plot(eta, reshape(MI(l,1,:), [1, length(eta)]),'LineWidth', 2);% tight-binding
hold on;
plot(eta, reshape(MI(l, length(theta),:), [1, length(eta)]), 'LineWidth', 2');% weak-binding
hold on;



y1 = 2+zeros(1,length(eta));
H1=area(eta,y1);
hold on
ide=eta>0;
H=area(eta(ide),y1(ide));
set(H(1),'FaceColor',[1 0.5 0]);
alpha(0.125)
set(gca,'FontSize',18)
legend({'tight-binding', 'weak-binding'}, 'Location', 'northwest');
xlabel('cross-talk level, $\eta$', 'Interpreter', 'latex');
ylabel('mutual information, $I$', 'Interpreter','latex');
title('$c=0$','Interpreter', 'latex', 'Fontsize', 15)
txt1 = 'no cross-talk';
txt2 = 'cross-talk';
t= text(3, 0.3, txt1, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
t= text(-8, 0.3, txt2, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
hold off;

print('MI_xtalk_zeroC' ,'-depsc');
print('MI_xtalk_zeroC' ,'-djpeg');




%% figure 4 -- worst case



hfig4 = figure(4);  
plot(eta, reshape(MI(l_minTB,1,:), [1, length(eta)]),'LineWidth', 2);% tight-binding
hold on;
plot(eta, reshape(MI(l_minWB, length(theta),:), [1, length(eta)]), 'LineWidth', 2');% weak-binding
hold on;


y1 = 2+zeros(1,length(eta));
H1=area(eta,y1);
hold on
ide=eta>0;
H=area(eta(ide),y1(ide));
set(H(1),'FaceColor',[1 0.5 0]);
alpha(0.125)
set(gca,'FontSize',18)
legend({'tight-binding', 'weak-binding'}, 'Location', 'northwest');
xlabel('cross-talk level, $\eta$', 'Interpreter', 'latex');
ylabel('mutual information, $I$', 'Interpreter','latex');
title(['$c_{TB}= $' num2str(Cr(l_minTB),2) ',  $c_{WB}= $' num2str(Cr(l_minWB),2)],'Interpreter', 'latex', 'Fontsize', 18)
txt1 = 'no cross-talk';
txt2 = 'cross-talk';
t= text(3, 1.5, txt1, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
t= text(-8, 1.5, txt2, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
hold off;


print('MI_xtalk_worst' ,'-depsc');
print('MI_xtalk_worst' ,'-djpeg');



%% figure 5-- x-tlak is most useful 
[v_maxTB, l_maxTB] = max(MI(:,1,1)-MI(:,1,length(eta)));
[v_maxWB, l_maxWB] = max(MI(:,length(theta),1)-MI(:,length(theta),length(eta)));

hfig5 = figure(5);  
plot(eta, reshape(MI(l_maxTB,1,:), [1, length(eta)]),'LineWidth', 2);% tight-binding
hold on;
plot(eta, reshape(MI(l_maxWB, length(theta),:), [1, length(eta)]), 'LineWidth', 2');% weak-binding
hold on;



y1 = 2+zeros(1,length(eta));
H1=area(eta,y1);
hold on
ide=eta>0;
H=area(eta(ide),y1(ide));
set(H(1),'FaceColor',[1 0.5 0]);
alpha(0.125)
set(gca,'FontSize',18)
legend({'tight-binding', 'weak-binding'}, 'Location', 'northwest');
xlabel('cross-talk level, $\eta$', 'Interpreter', 'latex');
ylabel('mutual information, $I$', 'Interpreter','latex');
title(['$c_{TB}= $' num2str(Cr(l_maxTB),2) ',  $c_{WB}= $' num2str(Cr(l_maxWB),2)],'Interpreter', 'latex', 'Fontsize', 18)
txt1 = 'no cross-talk';
txt2 = 'cross-talk';
t= text(3, 0.8, txt1, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
t= text(-8, 0.8, txt2, 'Interpreter', 'latex', 'Fontsize', 20);
t.VerticalAlignment = 'bottom';
hold off;

print('MI_xtalk_best' ,'-depsc');
print('MI_xtalk_best' ,'-djpeg');




%}
%% figure 6-- bar charts: positive, zero, and negative corr

figure(1); 

subplot(1,3,1) % c=0
y=[ MI(l,length(theta),length(eta)) MI(l,length(theta),1); ...
    MI(l,1,length(eta)) MI(l,1,1)];

c = categorical({'weak-binding','tight-binding'});
c = reordercats(c,{'weak-binding','tight-binding'});
h = bar(c,y);
set(gca,'FontSize',18)
ColOrd = get(gca,'ColorOrder');
h(2).FaceColor = ColOrd(5,:);
ylabel(' mutual information: \it I', 'Interpreter','tex');
legend('NO x-talk', 'x-talk','Location', 'northwest');
%title(['$c= $' num2str(Cr(l),4)], 'Interpreter','latex')
title([ 'correlation $c = 0$'], 'Interpreter','latex');


subplot(1,3,2) % c>0
y=[ mean(MI(lpos,length(theta),length(eta))) mean(MI(lpos,length(theta),1));...
    mean(MI(lpos,1,length(eta))) mean(MI(lpos,1,1))  ];

%c = categorical({'weak-binding','tight-binding'});
h = bar(c,y);
set(gca,'FontSize',18)
ColOrd = get(gca,'ColorOrder');
h(2).FaceColor = ColOrd(5,:);
%ylabel(' mutual info. $I$', 'Interpreter','latex');
%legend('x-talk', 'no x-talk','Location', 'best');
%title(['$c= $' num2str(Cr(lmax),4)], 'Interpreter','latex');
title([ 'correlation $c >0 $'], 'Interpreter','latex');

subplot(1,3,3)  % c<0
%lmin = 1271; % for this grid choice
y=[mean(MI(lneg,length(theta),length(eta))) mean(MI(lneg,length(theta),1)) ; ...
    mean(MI(lneg,1,length(eta)))  mean(MI(lneg,1,1))];

%c = categorical({'weak-binding','tight-binding'});
h = bar(c,y);
set(gca,'FontSize',18)
ColOrd = get(gca,'ColorOrder');
h(2).FaceColor = ColOrd(5,:);
%title(['$c= $' num2str(Cr(lmin),4)], 'Interpreter','latex');
title([ 'correlation $c <0$'], 'Interpreter','latex');



set(gcf, 'Position', [345, 565, 1305, 390]);
print('MI_xtalk_corr_bar' ,'-depsc');
print('MI_xtalk_corr_bar' ,'-djpeg');


% note the convention: corr_neg, cor_zero, cor_pos
ckt22WBnx = [mean(MI(lneg,length(theta),length(eta))) MI(l,length(theta),length(eta)) mean(MI(lpos,length(theta),length(eta)))];
ckt22WBwx = [mean(MI(lneg,length(theta),1)) MI(l,length(theta),1) mean(MI(lpos,length(theta),1))];
ckt22TBnx = [mean(MI(lneg,1,length(eta))) MI(l,1,length(eta)) mean(MI(lpos,1,length(eta)))];
ckt22TBwx = [mean(MI(lneg,1,1)) MI(l,1,1) mean(MI(lpos,1,1))];


%
ckt22WBnx_star = [MI(lcmin,length(theta),length(eta)) MI(l,length(theta),length(eta))  MI(lcmax,length(theta),length(eta))];
ckt22WBwx_star = [MI(lcmin,length(theta),1) MI(l,length(theta),1) MI(lcmax,length(theta),1)];
ckt22TBnx_star = [MI(lcmin,1,length(eta)) MI(l,1,length(eta)) MI(lcmax,1,length(eta))];
ckt22TBwx_star = [MI(lcmin,1,1) MI(l,1,1) MI(lcmax,1,1)];


newdata = [ckt22WBnx; ckt22WBwx; ckt22TBnx; ckt22TBwx;
    ckt22WBnx_star; ckt22WBwx_star; ckt22TBnx_star; ckt22TBwx_star];

csvwrite('two_by_twoDataSummary.csv', newdata);

fid = fopen('two_by_twoDataSummary.csv', 'wt');
fprintf(fid, '%s,%s,%s\n', 'corr_neg', 'corr_zero', 'corr_pos');
fprintf(fid, '%g,%g,%g\n', newdata.');   %transpose is important!
fclose(fid);

%{
%% figure 7-- bar charts: Pin for positive, zero, and negative corr
figure(7);

subplot(1,3,1)
c = categorical({'P(0,0)','P(1,0)','P(0,1)', 'P(1,1)'});
y = [Qall(l,:) 1-sum(Qall(l,:))];
barh(c,y,0.3);
set(gca,'FontSize',18)
title(['$c= $' num2str(Cr(l),4)], 'Interpreter','latex')

subplot(1,3,2)
c = categorical({'P(0,0)','P(1,0)','P(0,1)', 'P(1,1)'});
y = [Qall(lmax,:) 1-sum(Qall(lmax,:))];
barh(c,y,0.3);
set(gca,'FontSize',18)
title(['$c= $' num2str(Cr(lmax),4)], 'Interpreter','latex')

subplot(1,3,3)
c = categorical({'P(0,0)','P(1,0)','P(0,1)', 'P(1,1)'});
y = [Qall(lmin,:) 1-sum(Qall(lmin,:))];
barh(c,y,0.3);
set(gca,'FontSize',18)
title(['$c= $' num2str(Cr(lmin),4)], 'Interpreter','latex')


set(gcf, 'Position', [345, 565, 1305, 390]);

print('MI_xtalk_corr_pin' ,'-depsc');
print('MI_xtalk_corr_pin' ,'-djpeg');

%}
%% figure 8-- bar chat best worst
figure(2); 

y=[v_minWB,  v_maxWB;  ...
  v_minTB,  v_maxTB];

%c = categorical({'weak-binding','tight-binding'});
h = bar(c,y)
set(gca,'FontSize',18);
ColOrd = get(gca,'ColorOrder');
h(2).FaceColor = ColOrd(5,:);
ylabel('information gain: \it \Delta I', 'Interpreter','tex');
legend('least', 'most','Location', 'best');
set(gcf, 'Position', [680   564   412   414]);
print('MI_xtalk_best_worse_bar' ,'-depsc');
print('MI_xtalk_best_worse_bar' ,'-djpeg');






figure(3)
histogram(Cr,'Normalization','probability');
xlabel('input correlation, $c$', 'Interpreter', 'latex');
ylabel('$P(c)$', 'Interpreter','latex');
set(gca,'FontSize',28)
print('prob_of_corr' ,'-depsc');
print('prob_of_corr' ,'-djpeg');



save('xtalk_general_data');


   

%% function specification

function x = f(y)
x = exp(-y) /(1+exp(-y));
end

function x = g(y,z)

x = (exp(-y)+ exp(-z)) /(1+exp(-y)+exp(-z));

end
