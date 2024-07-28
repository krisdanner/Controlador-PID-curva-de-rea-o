%==============================================================
% Projeto de um controlador PID aplicado ao eixo z de um
% foguete pelo método da curva de reação em malha aberta 
%==============================================================
% Autor: Christian Danner Ramos de Carvalho
%==============================================================

close all
clear
clc

% Matrizes do sistema

A = [0 1;
    3.478e-09 9.103e-06];
B = [0.0388;
    0.9995];
C = [1 0];
D = 0;

eig(A)

M = ctrb(A,B);
rank(M)

N = obsv(A,C);
rank(N)

polos_desejados = [-1, -2];
K = place(A,B,polos_desejados);

% Sistema em malha fechada (alocação de pólos)

Amf = A-B*K;
Bmf = B;
Cmf = C;
Dmf = D;
plantaMF = ss(Amf, Bmf, Cmf, Dmf);

% Função de tranferência do sistema

[num,den] = ss2tf(Amf,Bmf,Cmf,Dmf);
G = tf(num,den);

% Tempo
dt = 0.05;
t = 0:dt:35;

% Entrada
u(1:length(t))=1;

% Saída

y = lsim(G,u,t);
%y = step(plantaMF,t2);
dy = diff(y)/dt; %Derivada

%Encontre o ponto de inflexão e sua derivada
% o ponto onde a inclinação da resposta ao degrau tem seu valor máximo (ponto de inflexão)

[m,p]=max(dy);
yp=y(p);
tp=t(p);
tm=0:7; 
ym=m*(tm-tp)+yp; %Equação da reta

% Plot da curva de reação e derivada no ponto de inflexão

figure
plot([t(1) t(end)],[0.5 0.5],'--k',[2.2 2.2],[0 0.5],'--k','linewidth',2);
grid
hold on
plot(t,y,tp,yp,'o',tm,ym,'-r','linewidth',3);
axis([0 35 0 0.6]); 
box off
ylabel('$$c(t)$$','FontSize',20,'Interpreter','latex')
xlabel('$$t$$','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',(20))

% Parametros obtidos pelo ponto de inflexão

K = 0.5;
L = 0.24;
T = 1.96;
G1 = tf([K],[T 1]);
G1.iodelay = L;
y1 = lsim(G1,u,t);

% Parâmetros obtidos pelo gráfico

L2 = 0.6;
T2 = 7/4;
G2 = tf([K],[T2 1]);
G2.iodelay = L2;
y2 = lsim(G2,u,t);

% Cálculo dos ganhos

Kp1 = T/(K*L);
Kp2 = 0.9*Kp1;
Kp3 = 1.2*Kp1;

Ki2 = Kp2/(L/0.3);
Ki3 = Kp3/(2*L);

Kd3 = Kp3*0.5*L;

% Plot comparação da resposta ao degrau por ponto de inflexão x ajuste
% pelo gráfico da função de transferência do processo:
% Gp(s) = (Ke^-Ls) / Ts + 1

figure
plot([t(1) t(end)],[0.5 0.5],'--k','linewidth',2);
hold on
plot(t,y,t,y1,'--',t,y2,'--','linewidth',3);
legend('Ganho K','Real','PID Inflexão', 'PID Gráfico');
%text(-1.5,2,'$$K$$','FontSize',20,'interpreter','latex')
axis([0 35 0 0.6]); 
box off
ylabel('$$c(t)$$','FontSize',20,'Interpreter','latex')
xlabel('$$t$$','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',(20))

% Plot PID por inflexão x PID por gráfico

Kc = (1.2*T)/(K*L);
Kc = Kc/1;  
ti = 2*L;
td = 0.5*L;

C = tf(Kc*[ti*td ti 1],[ti 0]);

Kc2 = (1.2*T2)/(K*L2);
Kc2 = Kc2/1;  
ti2 = 2*L2;
td2 = 0.5*L2;

C2 = tf(Kc2*[ti2*td2 ti2 1],[ti2 0]);

figure
step((C*G)/(1+C*G),(C2*G)/(1+C2*G))
legend('PID por Inflexão','PID por Gráfico')



