% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS 
%   + SEGUIMIENTO SE SETPOINT
clc; clear all; close all;

% N: tamaño del horizonte de predicción
% ejemplo: N = [2,3,4,5,10,20];
N = [10];

% Numero de iteraciones para PDIP
iterPDIP = 20;
% Linear solver: 'minres', 'cgrad', 'chol' o 'matlab'
linSol = "matlab";
% Número de iteraciones para minres
iterMINRES = 30;

% Comparar resultados con quadprog y LQR
compare = false;

% Guardar matrices de interés 
saveMat.MPC = true;
saveMat.PDIP = false;
saveMat.LS = false; % esta opción hace que sea mucho mas lento el código

% la función controlMpc necesita saber en que iteración va para
% guardar las matrices en el lugar que corresponde
saveMat.i = 0;

% Donde guardar las matrices que van a ser utilizadas como referencia
saveMat.outputMat = 'data/servoMats_N';

% Ts: Periodo de muestreo en segundos
Ts = 0.001;
% tsimu: Tiempo de simulación en segundos
tsimu = 10;
% t: Arreglo de tiempo
t=0:Ts:tsimu-Ts;

% Datos del servomotor en tiempo discreto
K=39.08/27.92;
tau=1/27.92; 
A = [exp(-Ts/tau),0; tau*(1-exp(-Ts/tau)),1];
B = [K*(1-exp(-Ts/tau));Ts*K+tau*K*(exp(-Ts/tau)-1)];
C = [0,1];
% x0: velocidad y posicion angular inicial
x0= [0;0];

% Datos de las restricciones
% umin, umax en Volts:
umin=-10;
umax=10;
% xmin, xmax: rad/s, rad
xmin=[-5;-2];
xmax=[5;2];

% Gamma y Omega se utilizan para MPC y LQR
Gamma=0.1;
Omega=C'*C;
[Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);

% yref: salida deseada del sistema
yref=[ones(1,length(t)/4),-ones(1,length(t)/4), ...
      ones(1,length(t)/4),-ones(1,length(t)/4)];


% Cálculo de la señal de control óptima y estado - control MPC
for n=N
    disp(['Procesando horizonte de tamaño: ', num2str(n)])
    
    % Reservar espacio memoria para guardar matrices ----------------------
    if saveMat.PDIP
        saveMat.h = zeros(n,1,length(t));
        saveMat.cx = zeros(6*n,1,length(t));
        saveMat.u_tilde = zeros(n,1,length(t));
    end
    if saveMat.LS
        saveMat.Ak = [];
        saveMat.bk = [];
        saveMat.zk = [];
    end
    
    
    % Control MPC ---------------------------------------------------------
    % Reservar espacio de memoria
    xPDIP = zeros(2,length(t));
    uPDIP = zeros(1,length(t));
    % Estado inicial
    xPDIP(:,1) = x0;

    % Preparar MPC
    [Acal,Ocal,Omg,H,Mx] = setupMPC(A,B,OmegaN,Omega,Gamma,n);

    for i=1:length(t)
        saveMat.i = i;
        [uPDIP(i),saveMat]=controlMPC(xPDIP(:,i),yref(i),xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Omg,H,Mx,n,'pdip',linSol,iterPDIP,iterMINRES, saveMat);
        % Cálculo del siguiente estado
        xPDIP(:,i+1)=A*xPDIP(:,i)+B*uPDIP(i);
    end
    xPDIP = xPDIP(:,1:end-1);

    
    % Guardar matrices de interés -----------------------------------------
    if saveMat.MPC || saveMat.PDIP || saveMat.LS
        % Si es que no existe el .mat, se crea
        saveMat.outputMat = saveMat.outputMat + string(n) + '.mat';
        if exist(saveMat.outputMat, 'file')
            save(saveMat.outputMat,'n','-append');
        else
            save(saveMat.outputMat,'n', '-v7.3');
        end
        matObj = matfile(saveMat.outputMat,'Writable',true);
    end
    if saveMat.LS
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES; 
        matObj.Ak = saveMat.Ak;
        matObj.bk = saveMat.bk;
        matObj.zk = saveMat.zk;
        writeLSSamples(n, saveMat.outputMat);
    end
    if saveMat.PDIP
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES; 
        matObj.H = H;
        matObj.Mx = Mx;
        matObj.h = saveMat.h;
        matObj.cx = saveMat.cx;
        matObj.u_tilde = saveMat.u_tilde;
        writePDIPSamples(n, saveMat.outputMat);
    end
    if saveMat.MPC
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES; 
        matObj.umin = umin;
        matObj.umax = umax;
        matObj.xmin = xmin;
        matObj.xmax = xmax;
        matObj.Acal = Acal;
        matObj.A = A;
        matObj.B = B;
        L=[eye(size(A,1))-A,-B;C,zeros(1,1)]; 
        L_inv = inv(L);
        matObj.L_invLast = L_inv(:,3);
        matObj.AcalOmgOcal = Acal'*Omg*Ocal;
        matObj.Acal = Acal;
        matObj.yref = yref;
        matObj.H = H;
        matObj.Mx = Mx;
        matObj.u = uPDIP;
        matObj.x = xPDIP;
        writeMPCSamples(n, saveMat.outputMat);
    end

    
    % Comparar resultados con Quadprog y LQR ------------------------------
    if compare
        % Reservar espacio de memoria
        xQuadProg = zeros(2,length(t));
        uQuadProg = zeros(length(t),1);

        % Estado inicial
        xQuadProg(:,1) = x0;

        for i=1:length(t)
            [uQuadProg(i),saveMat]=controlMPC(xQuadProg(:,i),yref(i),xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Omg,H,Mx,n,'quadprog',linSol,iterPDIP,iterMINRES,saveMat);
            % Cálculo del siguiente estado
            xQuadProg(:,i+1)=A*xPDIP(:,i)+B*uQuadProg(i);
        end
        xQuadProg = xQuadProg(:,1:end-1);
        % Evolución del sistema ante la entrada óptima obtenida con LQR:
        xLQR(:,1)=x0;
        L=[eye(size(A,1))-A,-B;C,zeros(1,1)]; 
        for i=1:length(t)
            bl=[zeros(size(A,1),1);yref(1,i)];
            infy=L\bl; 
            xinfy=infy(1:size(A,1));         
            uinfy=infy(end);
            uLQR(1,i)=-Linf*(xLQR(:,i)-xinfy)+uinfy;
            xLQR(:,i+1) = A*xLQR(:,i) + B*uLQR(1,i);
        end 
        xLQR=xLQR(:,1:end-1);

        % Gráficas de los estados, salida y de la señal de control óptima
        % -----------Gráficas de las estimaciones del estado --------------
        figure; grid on; hold on
        plot(t,xPDIP(1,:)','-r','LineWidth', 3) 
        plot(t,xQuadProg(1,:)','-b','LineWidth', 2) 
        plot(t,xLQR(1,:),'--g','LineWidth', 1);
        legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
        title('Estado 1 - Velocidad angular N:'+string(n));
        box on; xlabel('t'); ylabel('wt'); 
        grid on; box on;  hold off;

        figure; grid on; hold on;
        plot(t,xPDIP(2,:)','-r','LineWidth', 3) 
        plot(t,xQuadProg(2,:)','-b','LineWidth', 2) 
        plot(t,xLQR(2,:),'--g','LineWidth', 1);
        legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
        title('Estado 2 - Posición angular N:'+string(n)')
        ylim([-2,2]);
        box on; xlabel('t'); ylabel('theta'); 
        grid on; box on;  hold off;
        
        % -----------Grafico de la entrada --------------------------------
        figure; grid on; hold on;
        plot(t,uPDIP,'-r','LineWidth', 3);
        plot(t,uQuadProg,'-b','LineWidth', 2);
        plot(t,uLQR,'--g','LineWidth', 1); 
        legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
        title('Entrada u_t (Vm) N:'+string(n)');
        ylim([-10,10]);
        ylabel('u_t'); xlabel('t'); 
        grid on; box on;  hold off;
    end
end


    




