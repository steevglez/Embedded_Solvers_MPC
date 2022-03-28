%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS + SEGUIMIENTO SE SETPOINT (DENSE)
    clc; clear all; close all;
    %Datos del sistema
    K=39.08/27.92; 
    tau=1/27.92;
    T = 0.001;
    A = [exp(-T/tau),0; tau*(1-exp(-T/tau)),1];
    B = [K*(1-exp(-T/tau));T*K+tau*K*(exp(-T/tau)-1)];
    C = [0,1];
    x0= [0;0];
    
    %Datos de las restricciones
    N=5;                              Ns=10000;
    umin=-10;                         umax=10;
    xmin=[-5;-2];                     xmax=[5;2];
    %xmin=[-5;-0.8];                     xmax=[5;0.8];
    Gamma=0.1; Omega=C'*C;  
    [Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);
    yref=[ones(1,Ns/4),-ones(1,Ns/4),ones(1,Ns/4),-ones(1,Ns/4)];

    % Cálculo de la señal de control óptima y estado - control MPC
    [up,xp]=input_state_constraints_rt_aguero(x0,yref,xmin,xmax,umin,umax,A,B,C,OmegaN,Omega,Gamma,N,Ns,'pdip');
    
    % Evolución del sistema ante la entrada óptima obtenida con LQR:
    xr(:,1)=x0;
    L=[eye(size(A,1))-A,-B;C,zeros(1,1)]; 
    for i=1:1:Ns 
        bl=[zeros(size(A,1),1);yref(1,i)];
        infy=L\bl; 
        xinfy=infy(1:size(A,1));         
        uinfy=infy(end);
        ur(1,i)=-Linf*(xr(:,i)-xinfy)+uinfy;
        xr(:,i+1) = A*xr(:,i) + B*ur(1,i);
    end 
    xr=xr(:,1:end-1);

    % Gráficas de los estados, salida y de la señal de control óptima
    t=1:1:Ns;
    % -----------Gráficas de las estimaciones del estado ----------------
    figure; grid on; hold on;
    plot(t,xp(1,:)','-r','LineWidth', 3) 
    plot(t,xr(1,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','LQR');  
    title('Estado 1 - Velocidad angular');
    box on; xlabel('t'); ylabel('wt'); 
    grid on; box on;  hold off;
    
    figure; grid on; hold on;
    plot(t,xp(2,:)','-r','LineWidth', 3) 
    plot(t,xr(2,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','LQR');  
    title('Estado 2 - Posición angular')
    ylim([-2,2]);
    box on; xlabel('t'); ylabel('theta'); 
    grid on; box on;  hold off;
    %---------------------------------Entrada------------------------------
    
    figure; grid on; hold on;
    plot(t,up,'-r','LineWidth', 3);
    plot(t,ur,'--g','LineWidth', 1); 
    legend('MPC_{pdip}','LQR');  
    title('Entrada u_t (Vm)');
    ylim([-10,10]);
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  hold off;
    
function [u_opt,xt]=input_state_constraints_rt_aguero(x0,yref,xmin,xmax,umin,umax,A,B,C,OmegaN,Omega,Gamma,N,Ns,qpsol)
    
    options =  optimset('Display','off');
    [Acal,Ocal] = ABcal(A,B,N);
    Gmm = mdiag(Gamma,N);                  
    Omg = mdiag(Omega,OmegaN,N);
    % Matriz H del funcional de costo VN(x_0,vec{u})
    H=2*Gmm+2*Ocal'*Omg*Ocal;              
    H=(H+H')/2;
    % Matriz M de las restricciones lineales Mu<=c
    M=[Ocal;-Ocal];  
    % Convertir restricciones caja en Desigualdad
    Mx=[M;eye(N);-eye(N)];
    % Nueva x0 considerando que se desea seguir una referencia
    xt(:,1)=x0;
  
    for i=1:1:Ns
        [xinfy,uinfy] = Infinity(A,B,C,yref(1,i));
        %xinfy=0; uinfy=0;
        % Vectores a y b de las restricciones a<=u<=b
        a=ones(N,1)*(umin-uinfy);                   
        b=ones(N,1)*(umax-uinfy);
        x0nau=x0-xinfy;
        % Vector f del funcional de costo VN(x_0,vec{u})
        f=(2*x0nau'*Acal'*Omg*Ocal)';
        % Vector c de las restricciones lineales Mu<=c
        [Xmin,Xmax] = xboxes(N,xmin,xmax,xinfy);
        c=[(Xmax-Acal*x0nau);-(Xmin-Acal*x0nau)];
        % Convertir restricciones caja en Desigualdad
        cx=[c;b;-a];
        % Solución del problema de optimización u^opt
        if strcmp(qpsol,'pdip')
            unau=pdip(H,f,Mx,cx,30,1e-9,30,'lschol');
        end
        if strcmp(qpsol,'quadprog')
            unau=quadprog(H,f,Mx,cx,[],[],[],[],[],options);
        end
        
        %Recolectar datos para pruebas en FPGA
        ft(:,:,i)=f;
        binq(:,:,i)=cx;
        sqp(:,:,i)=unau;
        %------------------------------------
        ut=unau(1)+uinfy;
        u_opt(i,1)=ut;
        % Cálculo del siguiente estado inicial para optmizar VN(x_0,vec{u})
        xt(:,i+1)=A*xt(:,i)+B*ut;
        x0=xt(:,i+1); 
    end
    
    save(['H_',int2str(N)],'H');
    save(['ft_',int2str(N)],'ft');
    save(['Ginq_',int2str(N)],'Mx');
    save(['binq_',int2str(N)],'binq');
    save(['sqp_',int2str(N)],'sqp');
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

function [Xmin,Xmax] = xboxes(N,min,max,xinfy)
    Xmin=[]; Xmax=[]; 
    for i=1:1:N
        Xmin=[Xmin;min-xinfy];
        Xmax=[Xmax;max-xinfy];
    end
end

function [xinfy,uinfy] = Infinity(A,B,C,yref)
    L=[eye(size(A,1))-A,-B;C,zeros(1,1)];              
    bl=[zeros(size(A,1),1);yref];
    infy=L\bl; 
    xinfy=infy(1:size(A,1));         
    uinfy=infy(end);
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
