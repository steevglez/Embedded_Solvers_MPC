%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS (DENSE)
    clc; clear all; close all;
    %Datos del sistema
    A=[1,0,0,0,0,0,0.5,0,0,0,0,0;...
       0,1,0,0,0,0,0,0.5,0,0,0,0;...
       0,0,1,0,0,0,0,0,0.5,0,0,0;...
       0,0,0,1,0,0,0,0,0,0.5,0,0;...
       0,0,0,0,1,0,0,0,0,0,0.5,0;...
       0,0,0,0,0,1,0,0,0,0,0,0.5;...
       -1,0.5,0,0,0,0,1,0,0,0,0,0;...
       0.5,-1,0.5,0,0,0,0,1,0,0,0,0;...
       0,0.5,-1,0.5,0,0,0,0,1,0,0,0;...
       0,0,0.5,-1,0.5,0,0,0,0,1,0,0;...
       0,0,0,0.5,-1,0.5,0,0,0,0,1,0;...
       0,0,0,0,0.5,-1,0.5,0,0,0,0,1];
   B=[0,0,0;...
      0,0,0;...
      0,0,0;...
      0,0,0;...
      0,0,0;...
      0,0,0;...
      0.5,0,0;...
      -0.5,0,0;...
      0,0.5,0;...
      0,0,0.5;...
      0,-0.5,0;...
      0,0,-0.5];
   
    x0=0.1*ones(12,1);

    %Datos de las restricciones
    N=10;                             Ns=150;            Omega=eye(12);                       
    umin=-0.5*ones(3,1);              umax=0.5*ones(3,1);
    xmin=-4*ones(12,1);               xmax=4*ones(12,1);
    Nxmin=xmin;                       Nxmax=xmax;
    Gamma=eye(3);                     OmegaN=Omega;
    
    % Cálculo de la señal de control óptima y estado - control MPC
    [Acal, Ocal, ut,xt]=input_state_constraints_mpc(x0,xmin,xmax,Nxmin,Nxmax,umin,umax,A,B,OmegaN,Omega,Gamma,N,Ns,'pdip');
    
    % Gráficas de los estados, salida y de la señal de control óptima
    t=1:1:Ns;
    % -----------Gráficas de las estimaciones del estado ----------------
    n=size(A,1); figure; 
    for i=1:2:n
        grid on; hold on;
        plot(t,xt(i,:)','LineWidth', 2) 
        box on; xlabel('t'); ylabel('x'); 
        grid on; box on;  hold off;
    end
    %---------------------------------Entrada------------------------------
    m=size(B,2); figure; 
    for i=1:1:m
        grid on; hold on;
        plot(t,ut(i,:),'LineWidth', 2);
        ylabel('u_t'); xlabel('t'); 
        grid on; box on;  hold off;
    end
    
function [Acal, Ocal, u_opt,xt]=input_state_constraints_mpc(x0,xmin,xmax,Nxmin,Nxmax,umin,umax,A,B,OmegaN,Omega,Gamma,N,Ns,qpsol)
    
    [Acal,AN,Ocal,CN,Xmin,Xmax] = ABcal(A,B,N,xmin,xmax);
    Gmm = mdiag(Gamma,N);                  
    Omg = mdiag(Omega,OmegaN,N);
    % Matriz H del funcional de costo VN(x_0,vec{u})
    H=2*Gmm+2*Ocal'*Omg*Ocal;              
    H=(H+H')/2;
    % Matriz M de las restricciones lineales Mu<=c
    M=[CN;-CN;Ocal;-Ocal];
    % Vectores a y b de las restricciones a<=u<=b
    a=[]; b=[];
    for i=1:1:N
        a=[a;umin];
        b=[b;umax];
    end
    % Convertir restricciones caja en Desigualdad
    Mx=[M;eye(N*size(B,2));-eye(N*size(B,2))];
    xt(:,1)=x0;
    for i=1:1:Ns
        % Vector f del funcional de costo VN(x_0,vec{u})
        f=(2*x0'*Acal'*Omg*Ocal)';
        % Vector c de las restricciones lineales Mu<=c
        c=[(Nxmax-AN*x0);-(Nxmin-AN*x0);(Xmax-Acal*x0);-(Xmin-Acal*x0)];
        % Convertir restricciones caja en Desigualdad
        cx=[c;b;-a];
        % Solución del problema de optimización u^opt
        if strcmp(qpsol,'pdip')
            u=pdip(H,f,Mx,cx,30,1e-9,30,'cg_wip');
        end
        if strcmp(qpsol,'quadprog')
            u=quadprog(H,f,Mx,cx,[],[],[],[]);
        end
        u_opt(:,i)=u(1:size(B,2));
        % Cálculo del siguiente estado inicial para optmizar VN(x_0,vec{u})
        xt(:,i+1)=A*xt(:,i)+B*u(1:size(B,2));
        x0=xt(:,i+1);
    end
    xt=xt(:,1:end-1);
end
% Función para generar las matrices de Acal*x_0 + Ocal*vec{u}
function [Acal,AN,Ocal,CN,Xmin,Xmax] = ABcal(A,B,N,min,max)
    Acal=[]; Ocal=[]; Xmin=[]; Xmax=[];
    for i=1:1:N
        AB=[];
        Acal=[Acal;A^i];
        for j=1:1:i; AB=[AB,A^(i-j)*B]; end
        Fl=[AB,zeros(size(AB,1),(N-i)*size(B,2))];
        Ocal=[Ocal;Fl];
        if i==N; AN=A^i; CN=Fl;  end
        Xmin=[Xmin;min];
        Xmax=[Xmax;max];
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

