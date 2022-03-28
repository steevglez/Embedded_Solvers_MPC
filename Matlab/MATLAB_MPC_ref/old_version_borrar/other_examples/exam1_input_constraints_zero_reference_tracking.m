%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA (DENSE)
    clc; clear all; close all;
    %Datos del sistema 
    ts=0.004; T=0.06; K=0.15;
    A=[1,ts;0,1-ts/T];
    B=[0;ts*(K/T)];
    x0= [0.8;-0.4];

    %Datos de las restricciones
    N=40;                             Ns=100;
    Omega=[1.1e4,0;0,2.9e1];          umin=-100;          umax=100;
    Gamma=2.4e-1;
    [Linf,~,~]=dlqr(A,B,Omega,Gamma);
    OmegaN = Omega;

    % Cálculo de la señal de control óptima y estado - control MPC
    [Acal, Ocal, ut,xt]=input_state_constraints_mpc(x0,umin,umax,A,B,OmegaN,Omega,Gamma,N,Ns,'pdip');
    
    % Evolución del sistema ante la entrada óptima obtenida con LQR:
    xr(:,1)=x0;
    for i=1:1:Ns   
        ur(i,1)=-Linf*xr(:,i);
        xr(:,i+1) = A*xr(:,i) + B*ur(i,1);
    end 
    xr=xr(:,1:end-1);
    
    % Gráficas de los estados, salida y de la señal de control óptima
    t=0.005:0.005:0.5;
    % -----------Gráficas de las estimaciones del estado ----------------
    figure; grid on; hold on;
    plot(t,xt(1,:),'LineWidth', 2) 
    plot(t,xt(2,:),'LineWidth', 2) 
    box on; xlabel('t'); ylabel('x'); 
    box on;  hold off;
 
    
    %---------------------------------Entrada------------------------------
    figure; 
    %plot(t,ut(i,:),'LineWidth', 2);
    plot(t,ut,'LineWidth', 2);
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  
 
    
    
function [Acal, Ocal, u_opt,xt]=input_state_constraints_mpc(x0,umin,umax,A,B,OmegaN,Omega,Gamma,N,Ns,qpsol)
    options =  optimset('Display','off');
    [Acal,Ocal] = ABcal(A,B,N);
    Gmm = mdiag(Gamma,N);                  
    Omg = mdiag(Omega,OmegaN,N);
    % Matriz H del funcional de costo VN(x_0,vec{u})
    H=2*Gmm+2*Ocal'*Omg*Ocal;              
    H=(H+H')/2;
    % Vectores a y b de las restricciones a<=u<=b
    a=ones(N,1)*umin;                    
    b=ones(N,1)*umax;
    % Convertir restricciones caja en Desigualdad
    Mx=[eye(N);-eye(N)];
    xt(:,1)=x0;
    for i=1:1:Ns
        % Vector f del funcional de costo VN(x_0,vec{u})
        f=(2*x0'*Acal'*Omg*Ocal)';
        % Convertir restricciones caja en Desigualdad
        cx=[b;-a];
        % Solución del problema de optimización u^opt
        if strcmp(qpsol,'pdip')
            u=pdip(H,f,Mx,cx,20,1e-9,30,'matlab');
        end
        if strcmp(qpsol,'quadprog')
            u=quadprog(H,f,Mx,cx,[],[],[],[],[],options);
        end
        u_opt(i,1)=u(1);
        % Cálculo del siguiente estado inicial para optmizar VN(x_0,vec{u})
        xt(:,i+1)=A*xt(:,i)+B*u(1);
        x0=xt(:,i+1);
    end
    xt=xt(:,1:end-1);
end
% Función para generar las matrices de Acal*x_0 + Ocal*vec{u}
function [Acal,Ocal] = ABcal(A,B,N)
    Acal=[]; Ocal=[]; 
    for i=1:1:N
        AB=[];
        Acal=[Acal;A^i];
        for j=1:1:i; AB=[AB,A^(i-j)*B]; end
        Fl=[AB,zeros(size(AB,1),N-i)];
        Ocal=[Ocal;Fl];
    end
end
% Función para generar matrices diagonales por bloques
function [M] = mdiag(varargin) 
    if nargin==2
        A=varargin{1}; N=varargin{2}; M=A;    
        for i=1:1:N-1; M=blkdiag(M,A); end
    end
    if nargin==3
        A=varargin{1}; B=varargin{2}; N=varargin{3};  M=A;
        for i=1:1:N-2; M=blkdiag(M,A); end
        M=blkdiag(M,B);
    end
end

