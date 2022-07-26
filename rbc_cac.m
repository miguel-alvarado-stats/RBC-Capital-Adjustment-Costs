% rbc_cac

% Miguel Alvarado

clear global;
clear all;
clc;
% cd D:\matlab_my_files\malvarado_programs\ejercicio2_a;
format compact;
format short;

global p_beta p_delta p_alfa p_sigma p_v p_chi p_rho ...
       c_ss k_ss z_ss n_ss y_ss i_ss w_ss R_ss q_ss lambda_ss;

% parametros

p_beta  = 0.9825;
p_delta = 0.025;
p_alfa  = 1/3;
p_sigma = 2;
p_v     = 1;
p_chi   = 1;
p_rho   = 0.9;

% Calculo Estado EStacionario (ver archivo ss_for_seq_rbc_cac.m)

run ss_for_seq_rbc_cac.m

% Para el sistema log-linealizado colapsado (ver archivo pdf, Sistema 2), 
% denotamos los siguientes parámetros

phi_1 = ((p_sigma*(1-p_alfa))/(p_v+p_alfa));
phi_2 = ((p_v*(p_alfa-1))/(p_v+p_alfa));
phi_3 = ((p_alfa*(p_v+1))/(p_v+p_alfa));
phi_4 = ((p_v+1)/(p_v+p_alfa));

m_1 = y_ss/i_ss;
m_2 = c_ss/i_ss;

a_1 = -(p_sigma+((p_chi*p_delta)*((phi_1*m_1)+m_2)));
a_2 = ((phi_3*m_1)-1)*p_chi*p_delta;
a_3 = p_chi*p_delta*(phi_4*m_1);

b_1 = -(p_sigma+(p_beta*R_ss*phi_1)+((p_beta*p_chi*p_delta)*((phi_1*m_1)+m_2)));
b_2 = (p_beta*R_ss*phi_2)+(p_beta*((phi_3*m_1)-1)*p_chi*p_delta);
b_3 = (p_beta*R_ss*phi_4)+(p_beta*p_chi*p_delta*(phi_4*m_1));

d_1 = -p_delta*((phi_1*m_1)+m_2);
d_2 = (p_delta*(phi_3*m_1))+(1-p_delta);
d_3 = p_delta*(phi_4*m_1);

% modelo matricial;

A1 = [a_1 a_2 a_3;
      d_1 d_2 d_3;
      0 0 p_rho];
      
A2 = [b_1 b_2 b_3;
      0 1 0;
      0 0 1];
  
A3 = [0 b_1 b_2 b_3;
      0 0 0 0;
      -1 0 0 0];
  
% formal estructural  
 
A1_inv = inv(A1);

A = A1_inv*A2;

B = A1_inv*A3;

[Q,F] = jordan(A); %F matriz de autovalores de A, Q matriz de autovectores de A

[U,T] = schur(A);  %T matriz de autovalores de A, U matriz de autovectores de A

% Verificamos las condiciones de Blanchard & Kahn (estabilidad del Sistema 2)

variables = 3;
vcontrol = 1; %variablles libres {c_(t)}

root_g_1=0; % contador
root_l_1=0; % contador

for i=1:variables
if abs(F(i,i))>=1
    root_g_1=root_g_1 + 1;
else
    root_l_1=root_l_1 + 1;
end
end;

if root_l_1==vcontrol
    display('SISTEMA ESTABLE'), F
else
    display('SISTEMA INESTABLE'), break
end

% Funciones de Politica (FP) & IRF

QQ=Q^-1;

q=zeros(root_l_1,variables); % vector donde recogeré los coeficientes para armar la FP-IRF

j=1; % contador

for i=1:variables
    if abs(F(i,i))<1
        q(j,:)=QQ(i,:);
        j=j+1;
    end
end

% vectores donde recogeré la trayectoria de las variables del modelo

z_tray = zeros(50,1);
c_tray = zeros(50,1);
k_tray = zeros(51,1);

% Dinámica del shock de 1% en el AR(1) del modelo

z_tray(1,1) = 1; % Shock en t = 1
k_tray(1,1) = 0; % Cero al ser variable predeterminada (estado) en t = 1.
c_tray(1,1) = -(q(1,2)/q(1,1))*k_tray(1,1) - (q(1,3)/q(1,1))*z_tray(1,1); % Ecuación (14) en t = 1

for i=2:50
    z_tray(i,1) = p_rho*z_tray(i-1,1); % Ecuación (16)
    k_tray(i,1) = (d_1*c_tray(i-1,1)) + (d_2*k_tray(i-1,1)) + (d_3*z_tray(i-1,1)); % Del sistema log-linealizado: % Ecuación (15)
    c_tray(i,1) = -(q(1,2)/q(1,1))*k_tray(i,1) - (q(1,3)/q(1,1))*z_tray(i,1); % Ecuación (14)
end

% calculamos un periodo adicional para el stock de capital.

k_tray(51,1) = (d_1*c_tray(50,1)) + (d_2*k_tray(50,1)) + (d_3*z_tray(50,1));

% calculamos el resto de las trayectorias:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores donde recogeré la trayectoria de las variables restantes del modelo

lambda_tray = zeros(50,1);
n_tray  = zeros(50,1);
i_tray  = zeros(50,1);
y_tray  = zeros(50,1);
w_tray = zeros(50,1);
R_tray = zeros(50,1);
q_tray  = zeros(50,1);

for i=1:50
    lambda_tray(i,1) = (-1)*p_sigma*c_tray(i,1); % Ecuacion (1)
    n_tray(i,1) = (-p_sigma/(p_v+p_alfa))*c_tray(i,1) + (1/(p_v+p_alfa))*z_tray(i,1) + (p_alfa/(p_v+p_alfa))*k_tray(i,1); % Ecuacion (17)
    w_tray(i,1) = z_tray(i,1) + p_alfa*k_tray(i,1) - p_alfa*n_tray(i,1); % Ecuacion (5)
    R_tray(i,1) = z_tray(i,1) + (p_alfa - 1)*k_tray(i,1) + (1 - p_alfa)*n_tray(i,1); % Ecuacion (6)
    y_tray(i,1) = z_tray(i,1) + p_alfa*k_tray(i,1) + (1 - p_alfa)*n_tray(i,1); % Ecuacion (7)
    i_tray(i,1) = (m_1*y_tray(i,1)) - (m_2*c_tray(i,1)); % Ecuacion (8)
    q_tray(i,1) = (p_chi*p_delta)*(i_tray(i,1)-k_tray(i,1)); % Ecuacion 3
end


% Ajusto el lag que realiza dynare en el stock de capital

k_tray = [k_tray(2:51,1)];

% Variables para los gráficos.

time = [0:1:49];
l_cero = zeros(50,1);

% Grafico de las IRF's

% Trayectoria tecnología: z
subplot(3,4,1), plot(time,z_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(z_tray)+0.05]), title('Tecnología: z')

% Trayectoria Consumo: c
subplot(3,4,2), plot(time,c_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(c_tray)+0.05]), title('Consumo: c')

% Trayectoria Stock de Capital: k
subplot(3,4,3), plot(time,k_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(k_tray)+0.05]), title('Stock de Capital: k')

% Trayectoria Trabajo: n
subplot(3,4,4), plot(time,n_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(n_tray)-0.05 max(n_tray)+0.05]), title('Trabajo: n')

% Trayectoria Inversión: i
subplot(3,4,5), plot(time,i_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(i_tray)-0.05 max(i_tray)+0.05]), title('Inversión: i')

% Trayectoria Producto: y
subplot(3,4,6), plot(time,y_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(y_tray)-0.05 max(y_tray)+0.05]), title('Producto: y')

% Trayectoria Salarios: w
subplot(3,4,7), plot(time,w_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(w_tray)+0.05]), title('Salarios: w')

% Trayectoria Costo del Capital: R
subplot(3,4,8), plot(time,R_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(R_tray)-0.05 max(R_tray)+0.05]), title('Costo del Capital: R')

% Trayectoria Lambda: Lambda
subplot(3,4,9), plot(time,lambda_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(lambda_tray)-0.1 0.05]), title('Lambda: Lambda')

% Trayectoria q de Tobin: q
subplot(3,4,10), plot(time,q_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(q_tray)-0.05 max(q_tray)+0.05]), title('q de Tobin: q')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ESTA ES OTRA SOLUCION EQUIVALENTE %%%%%%%

A1X = [-p_sigma -p_chi 0;
      d_1 d_2 d_3;
      0 0 p_rho];
  
A2X = [b_1 (b_2-p_chi) b_3;
       0 1 0;
       0 0 1];

A1X_inv = inv(A1X);   

AX = A1X_inv*A2X;

[QX,FX] = jordan(AX);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%