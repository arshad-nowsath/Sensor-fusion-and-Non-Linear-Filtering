clc;
clear all;
close all;
A=1;
Q=1.5;
R=2.5;
x_0=2;
P_0=6;
N=20;
H=1;
X = genLinearStateSequence(x_0, P_0, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);
%correlation
figure(1);
plot(X(2:N+1),Y,'*');
axis equal
%measurement flux
figure(2)
plot([1:N],X(2:N+1),'-',[1:N],Y,'*');
hold on, grid on;
xlabel 'time step', ylabel 'value'
legend('true','meaesurement')
%% 1b
[Xfiltered, P,xp,Pp,l] = kalmanFilterextract(Y, x_0, P_0, A, Q, H, R);
Xfilteredp3sig=[Xfiltered] + 3*sqrt([P(:)']);
Xfilteredm3sig=[Xfiltered] - 3*sqrt([P(:)']);
figure(3)
plot([1:N],Xfiltered,'-',[1:N],X(2:N+1),'-',[1:N],Y,'*',[1:N],Xfilteredp3sig,'-',[1:N],Xfilteredm3sig,'-');
hold on, grid on;
xlabel 'time step', ylabel 'value'
legend('filtered est','true','meaesurement','+3sigma','-3sigma')

timeframe = [5,10,15,20];
for i=1:length(timeframe)
    xmeanest  = Xfiltered(timeframe(i));
    xvar = P(timeframe(i));
    xnormpdf = xmeanest-(4*xvar):0.01:xmeanest+(4*xvar);
    ynormpdf = normpdf(xnormpdf,xmeanest,sqrtm(xvar));
    figure()
    plot(xnormpdf,ynormpdf, 'LineWidth',2);
    hold on
    xtrue = X(timeframe(i)+1);
    plot([xtrue,xtrue], [0,max(ynormpdf)],'--', 'LineWidth',3);   
    ylim([0,max(ynormpdf)]);
    title(['posterior distribution for t']);
    legend('PDF','true')
end
%% 1C
timeframe1c=3
    
    xmeanest_prior  = Xfiltered(timeframe1c-1);
    xvar_prior = P(timeframe1c-1);
    xnormpdf_prior = xmeanest_prior-(4*xvar_prior):0.01:xmeanest_prior+(4*xvar_prior);
    ynormpdf_prior = normpdf(xnormpdf_prior,xmeanest_prior,sqrtm(xvar_prior));

    xmeanest_prediict  = xp(timeframe1c);
    xvar_prediict = Pp(timeframe1c);
    xnormpdf_prediict = xmeanest_prediict-(4*xvar_prediict):0.01:xmeanest_prediict+(4*xvar_prediict);
    ynormpdf_prediict = normpdf(xnormpdf_prediict,xmeanest_prediict,sqrtm(xvar_prediict));
    
    xmeanest_posterior  = Xfiltered(timeframe1c);
    xvar_posterior = P(timeframe1c);
    xnormpdf_posterior = xmeanest_posterior-(4*xvar_posterior):0.01:xmeanest_posterior+(4*xvar_posterior);
    ynormpdf_posterior = normpdf(xnormpdf_posterior,xmeanest_posterior,sqrtm(xvar_posterior));
    
    figure()
    plot(xnormpdf_prior,ynormpdf_prior,'-', xnormpdf_prediict,ynormpdf_prediict,'-',xnormpdf_posterior,ynormpdf_posterior ,'-')
    hold on
    xtrue1 = X(timeframe1c);
    plot([xtrue1,xtrue1], [0,max(ynormpdf_prior)*1.2],'--', 'LineWidth',3);   
    hold on
    measure = Y(timeframe1c);
    plot([measure,measure], [0,max(ynormpdf_prior)*1.2],'--', 'LineWidth',3);  
    legend('prior','prediction','posterior','true','measurement')
    
   %% 1) D)
close all;
clear all;
A=1;
Q=1.5;
R=2.5;
x_0=2;
P_0=6;

H=1;
N = 5000;
X = genLinearStateSequence(x_0,P_0,A,Q,N);
Y = genLinearMeasurementSequence(X, H, R);


[Xfiltered, P,xp,Pp,l] = kalmanFilterextract(Y, x_0, P_0, A, Q, H, R);

Xfilteredm=mean(Xfiltered);
vm=mean([l.innov])
% 
  xnormpdf = 0-(4*P(:,:,end)):0.01:0+(4*P(:,:,end));
  ynormpdf = normpdf(xnormpdf,0,sqrtm(P(:,:,end)))
 figure();
 hold on, grid on;
 histogram( (Xfiltered-X(:,2:end)), 100 ,'Normalization','pdf');
 hold on,
 plot(xnormpdf,ynormpdf, 'LineWidth',2 );
 xlabel('$x_k-\hat{x}_{k|k}$','Interpreter','Latex')
 ylabel('normalized frequency','Interpreter','Latex')
 title 'Histogram of normalized estimation error'
 legend('histogram','Probablity density function')

 figure();
 hold on, grid on;A
 autocorr([l.innov])
 xlabel 'Lag', ylabel 'autocorrelation', title ' autocorrelation of innovation'
%% 1 f
clc;
close all;
clear all;
A=1;
Q=1.5;
R=2.5;
x_0=2;
P_0=6;
x_w=10;
H=1;
N = 20;
X = genLinearStateSequence(x_0,P_0,A,Q,N);
Y = genLinearMeasurementSequence(X, H, R);


[Xfiltered, P,xp,Pp,l] = kalmanFilterextract(Y, x_0, P_0, A, Q, H, R);


[Xfilteredw, Pw,xpw,Ppw,lw] = kalmanFilterextract(Y, x_w, P_0, A, Q, H, R);

figure()
plot([0:N],[x_0 Xfiltered],'-',[0:N],[x_0 X(2:N+1)],'-',[1:N],Y,'*',[0:N],[x_w Xfilteredw],'-');
hold on, grid on;
xlabel 'time step', ylabel 'value'
legend('filtered est','true','meaesurement','wrong prior filtered est')
%% 2A
clc;
clear all;
close all;
T=0.01;
A2=[1 T;0 1];
V2=[0;1];
H2=[1 0];
N2= 1;
muQ=[0;0];
varQ=[0 0;0 1.5];
R=2;
x_0=[1 ;3];
P_0=4*eye(2);
N=500;
X = genLinearStateSequence(x_0,P_0,A2,varQ,N);
Y = genLinearMeasurementSequence(X, H2, R);

figure()
plot([1:N],X(1,2:N+1),'-',[1:N],Y,'*');
hold on, grid on;
xlabel 'time step', ylabel 'value'
legend('true','meaesurement')

figure()
plot([1:N],X(2,2:N+1),'-');
grid on;
xlabel 'time step', ylabel 'value'
legend('Velocity')
%% 2B

[Xfiltered, P,xp,Pp,l] = kalmanFilterextract(Y, x_0, P_0, A2, varQ, H2, R);

Xfilteredp3sig=[Xfiltered] + 3*sqrt([squeeze(P(1,1,:))']);
Xfilteredm3sig=[Xfiltered] - 3*sqrt([squeeze(P(1,1,:))']);
figure()
plot([1:N],Xfiltered(1,:),'-',[1:N],X(1,2:N+1),'-',[1:N],Y,'*',[1:N],Xfilteredp3sig(1,:),'-',[1:N],Xfilteredm3sig(1,:),'-');
hold on, grid on;
xlabel 'time step', ylabel 'distance/position'
legend('filtered est','true','measurement','+3sigma','-3sigma')
figure()
plot([1:N],Xfiltered(2,:),'-',[1:N],X(2,2:N+1),'-',[1:N],Xfilteredp3sig(2,:),'-',[1:N],Xfilteredm3sig(2,:),'-');
hold on, grid on;
xlabel 'time step', ylabel 'Velocity'
legend('filtered est','true','+3sigma','-3sigma')
%% 2C
clc;
clear all;
close all;
T=0.01;
A2=[1 T;0 1];
V2=[0;1];
H2=[1 0];
N2= 1;
muQ=[0;0];
varQi=[0 0;0 1.5]
Q=[0.1 1 10 15];
R=2;
x_0=[1 ;3];
P_0=4*eye(2);
N=500;
X = genLinearStateSequence(x_0,P_0,A2,varQi,N);
Y = genLinearMeasurementSequence(X, H2, R);
for i=1:4
varQ=[0 0;0 Q(i)];

[G(i).Xfiltered, P,xp,Pp,l] = kalmanFilterextract(Y, x_0, P_0, A2, varQ, H2, R);
end
figure()
plot([1:N],G(1).Xfiltered(1,:),'-',[1:N],G(2).Xfiltered(1,:),'-',[1:N],G(3).Xfiltered(1,:),'-',[1:N],G(4).Xfiltered(1,:),'-');
hold on, grid on;
xlabel 'time step', ylabel 'distance/position'
legend('filtered est dist Q1','filtered est dist Q2','filtered est dist Q3','filtered est dist Q4')
figure()
plot([1:N],G(1).Xfiltered(2,:),'-',[1:N],G(2).Xfiltered(2,:),'-',[1:N],G(3).Xfiltered(2,:),'-',[1:N],G(4).Xfiltered(2,:),'-');
hold on, grid on;
xlabel 'time step', ylabel 'Velocity'
legend('filtered est dist Q1','filtered est dist Q2','filtered est dist Q3','filtered est dist Q4')
