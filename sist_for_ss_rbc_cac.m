% sist_for_ss_rbc_cac

% Miguel Alvarado

% Escribimos una funcion (sist_for_ss_rbc_cac), para nuestro sistema no 
% lineal deterministico F = 0, que determina el Estado Estacionario (EE) 
% de las variables de interes.
% El sistema viene de las CPO, tomamos las variables en estado estacionario
% esto es: x_(t) = x_(t+1) = x , for all t.
% Esta función acepta por input un vector X, donde X es un vector de variables 
% desconocidas (los EE) que luego son asignadas a una 
% respectiva variable del sistema (no importa el orden) y que genere (la 
% funcion) un vector F(X)como output.

% Además, añadimos como argumento de la función, un otro input P, donde P 
% es un vector de parámetros que se definen abajo de manera general
% (o bien, se definen abajo sus valores numericos).

function F=sist_for_ss_rbc_cac(X, P)

c      = X(1);
k      = X(2);
z      = X(3);
n      = X(4);
y      = X(5);
i      = X(6);
w      = X(7);
R      = X(8);
q      = X(9);
lambda = X(10);     % lambda

p_beta  = P(1);
p_delta = P(2);
p_alfa  = P(3);
p_sigma = P(4);
p_v     = P(5);
p_chi   = P(6);
p_rho   = P(7);

% (Ver archivo pdf, Sistema 1)

F = [(c^(-p_sigma)) - lambda;
     -(n)^(p_v) + lambda*w;
     q - 1;
     (q/p_beta) - R - (1-p_delta)*q;
     w - (1-p_alfa)*(k^p_alfa)*(n^-p_alfa);
     R - p_alfa*(k^(p_alfa-1))*(n^(1-p_alfa));
     y - (k^p_alfa)*(n^(1-p_alfa));
     y - c - i;
     i - p_delta*k;
     z - 1];
end
    
    