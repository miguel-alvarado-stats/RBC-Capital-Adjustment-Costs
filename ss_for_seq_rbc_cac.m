% ss_for_seq_rbc_cac

% Miguel Alvarado

% Para hallar el estado estacionario del sistema de ecuaciones definido
% en el m.file "sist_for_ss_rbc_cac.m" 

clear;
clc;

global p_beta p_delta p_alfa p_sigma p_v p_chi p_rho ...
       c_ss k_ss z_ss n_ss y_ss i_ss w_ss R_ss q_ss lambda_ss;
   
% Definimos el vector de parámetros P

P = [p_beta; p_delta; p_alfa; p_sigma; p_v; p_chi; p_rho];

% Estado Estacionario

opt = optimset('Display','iter');

% Tomaremos un valor inicial (guess) para todas las variables del modelo

X0=[0.5;
    0.5;
    0.5;
    0.5;
    0.5;
    0.5;
    0.5;
    0.5;
    0.5;
    0.5];

[solstar] = fsolve(@(X) sist_for_ss_rbc_cac(X,P),X0,opt);

% La solución hallada será, considerando la construcción de X en la función
% sist_for_ss_rbc_cac, se tiene:

c_ss = solstar(1,1);
k_ss = solstar(2,1);
z_ss = solstar(3,1);
n_ss = solstar(4,1);
y_ss = solstar(5,1);
i_ss = solstar(6,1);
w_ss = solstar(7,1);
R_ss = solstar(8,1);
q_ss = solstar(9,1);
lambda_ss = solstar(10,1);

display (c_ss);
display (k_ss);
display (z_ss);
display (n_ss);
display (y_ss);
display (i_ss);
display (w_ss);
display (R_ss);
display (q_ss);
display (lambda_ss);

