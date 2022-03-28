%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS + SEGUIMIENTO SE SETPOINT (SPARSE)
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
    N=5;                               Ns=10000;
    umin=-10;                          umax=10;
    xmin=[-5;-2];                      xmax=[5;2];
    %xmin=[-5;-0.8];                     xmax=[5;0.8];
    Gamma=0.1; Omega=C'*C;  
    [Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);
    yref=[ones(1,Ns/4),-ones(1,Ns/4),ones(1,Ns/4),-ones(1,Ns/4)];

    % Cálculo de la señal de control óptima y estado - control MPC
    [ut,xt]=input_state_constraints_rt_jerez(x0,yref,xmin,xmax,umin,umax,A,B,C,OmegaN,Omega,Gamma,N,Ns,'pdip');
    
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
    plot(t,xt(1,:)','-r','LineWidth', 2) 
    plot(t,xr(1,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','LQR');  
    title('Estado 1 - Velocidad angular');
    box on; xlabel('t'); ylabel('wt'); 
    grid on; box on;  hold off;
    
    figure; grid on; hold on;
    plot(t,xt(2,:)','-r','LineWidth', 2) 
    plot(t,xr(2,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','LQR');  
    title('Estado 2 - Posición angular')
    ylim([-2,2]);
    box on; xlabel('t'); ylabel('theta'); 
    grid on; box on;  hold off;
    %---------------------------------Entrada------------------------------
    
    figure; grid on; hold on;
    plot(t,ut,'-r','LineWidth', 2);
    plot(t,ur,'--g','LineWidth', 1); 
    legend('MPC_{pdip}','LQR');  
    title('Entrada u_t (Vm)');
    ylim([-10,10]);
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  hold off;
    
function [ut,xt]=input_state_constraints_rt_jerez(x0,yref,xmin,xmax,umin,umax,A,B,C,QN,Q,R,N,Ns,qpsol)   

    options =  optimset('Display','off');
% Datos y restricciones del sistema    
    [n,~]=size(A);  In=eye(n); 
    [~,m]=size(B);  Im=eye(m); 
    Znm=zeros(n,m); Zmn=zeros(m,n); 
    IN=eye(N); M=zeros(n,m); 
    J=[In;-In;Zmn;Zmn];       
    JN=[In;-In];             
    E=[Znm;Znm;Im;Im];
% Matrices para el problema de optimización cuadrática
    HT=kron(IN,[Q,M;M',R]);
    H = complete(HT,QN);
    H=(H+H')/2;
    h = zeros(size(H,1),1);
    F = mdiag(A,B,N); 
    GT=kron(IN,[J,E]);
    G = complete(GT,JN);
    Mx=[G;F;-F];
    xt(:,1)=x0;
    for i=1:1:Ns 
        [xinfy,uinfy] = Infinity(A,B,C,yref(1,i));
        d=[xmax-xinfy;-(xmin-xinfy);umax-uinfy;-(umin-uinfy)];
        dN=[xmax-xinfy;-(xmin-xinfy)];
        % Vector f de las restricciones de igualdad F*theta=f
        x0nau=x0-xinfy;
        f=[-x0nau;zeros(size(F,1)-n,1)];
        % Vector g de las restricciones de desigualdad G*theta<=g
        g=mvec(N,d,dN);
        % Convertir las restricciones de igualdad y desigualdad en
        % restricciones solo de desigualdad
        cx=[g;f;-f];
        % Solución del problema de optimización u^opt
        if strcmp(qpsol,'pdip')
            unau=pdip(H,h,Mx,cx,20,1e-9,50,'minres');
        end
        if strcmp(qpsol,'quadprog')
            unau=quadprog(H,h,Mx,cx,[],[],[],[],[],options);
        end
        
        %Recolectar datos para pruebas en FPGA
        ft(:,:,i)=h;
        binq(:,:,i)=cx;
        sqp(:,:,i)=unau;
        %------------------------------------
        
        up=unau(n+1:n+m)+uinfy;
        ut(:,i)=up;
        % Cálculo del siguiente estado inicial para optmizar VN(x_0,vec{u})
        xt(:,i+1)=A*xt(:,i)+B*up;
        yt(:,i)=C*xt(:,i+1);
        x0=xt(:,i+1);
    end
    xt=xt(:,1:end-1);
    
    save(['JH_',int2str(N)],'H');
    save(['Jft_',int2str(N)],'ft');
    save(['JGinq_',int2str(N)],'Mx');
    save(['Jbinq_',int2str(N)],'binq');
    save(['Jsqp_',int2str(N)],'sqp');
 
end

function [F] = mdiag(A,B,N)
    n=size(A,1);        m=size(B,2);         In=eye(n);
    f=[A,B,-In];        F=[];                col_F=(N+1)*n+N*m;      
    for i=1:1:N
        if i==1; F=[-In,zeros(n,col_F-n)]; end
        M=[zeros(n,(i-1)*(n+m)),f,zeros(n,col_F-((i-1)*(n+m)+(2*n+m)))]; F=[F;M];
    end
end

function [D] = complete(A,B)
    [fA,cA] = size(A);       [fB,cB] = size(B);
    Z1=zeros(fA,cB);          Z2=zeros(fB,cA);
    D=[A,Z1;Z2,B];
end

function [g]=mvec(N,d,dT)
    g=[];
    for jm=1:1:N+1
        if jm~=N+1; g=[g;d]; else; g=[g;dT]; end
    end
end

function [xinfy,uinfy] = Infinity(A,B,C,yref)
    L=[eye(size(A,1))-A,-B;C,zeros(1,1)];              
    bl=[zeros(size(A,1),1);yref];
    infy=L\bl; 
    xinfy=infy(1:size(A,1));         
    uinfy=infy(end);
end
