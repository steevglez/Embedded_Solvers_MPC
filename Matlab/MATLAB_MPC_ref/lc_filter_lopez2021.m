%%  %% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS - FILTRO LC 
    clc; clear all; close all;
    
    % Datos del sistema
    x0= [0;0;0;0];
    A =[0.5832    0.0367   -0.0570   -0.0036;
        -0.0367    0.5832    0.0036   -0.0570;
        11.3967    0.7167    0.5869    0.0369;
        -0.7167   11.3967   -0.0369    0.5869];
    B =[0.0571    0.0016;
       -0.0016    0.0571;
        0.4115    0.0170;
       -0.0170    0.4115];
    C=[0 0 1 0
       0 0 0 1];
       
    % Datos de las restricciones
    N=5;                              Ns=100;
    umin=[-420;-420];                 umax=[420;420];
    xmin=[-20;-20; -630; -630];       xmax=[10;20; 630; 630];
    Gamma=0.01*eye(2);                Omega=eye(4); 
    [Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);
    yref=[260*ones(1,Ns);zeros(1,Ns)];
    
    
    % Systems information
    G.sys.A=A;
    G.sys.B=B;
    G.sys.C=C;
    G.sys.x0=x0;
    G.sys.Ns=Ns;
    G.sys.yref=yref;
    
    % Cost function information
    G.cost.N=N;
    G.cost.OmegaN=OmegaN;
    G.cost.Omega=Omega;
    G.cost.Gamma=Gamma;
    
    % Constraints information
    G.cons.xmin=xmin;
    G.cons.xmax=xmax;
    G.cons.umin=umin;
    G.cons.umax=umax;
    
    % Linear solver information
    G.linsol.solver='minres';
    G.linsol.minres_iter=30;
    G.linsol.tol=1e-9;
    
    % QP solver information
    G.qpsol.solver='pdip';
    G.qpsol.pdip=20;
    
    % Cálculo de la señal de control óptima y estado - control MPC
    [ut,xp,~]=general_dense_mpc_instate_rt(G);
    
    G.qpsol.solver='quadprog';
    [up,xt,~]=general_dense_mpc_instate_rt(G);
    
    % Evolución del sistema ante la entrada óptima obtenida con LQR:
    L=[eye(size(A,1))-A,-B;C,zeros(2,2)];              
    
    xr(:,1)=x0;
    for i=1:1:Ns  
        bl=[zeros(size(A,1),1);yref(:,i)];
        infy=L\bl; 
        xinfy=infy(1:size(A,1));         
        uinfy=infy(size(A,1)+1:end);
        ur(:,i)=-Linf*(xr(:,i)-xinfy)+uinfy;
        xr(:,i+1) = A*xr(:,i) + B*ur(:,i);
        yr(:,i)=C*xr(:,i);
    end 
    xr=xr(:,1:end-1);
    
    % Gráficas de los estados, salida y de la señal de control óptima
    t=1:1:Ns;
    % -----------Gráficas de las estimaciones del estado ----------------
    figure; subplot(4,1,1); grid on; hold on;
    plot(t,xp(1,:)','-r','LineWidth', 3) 
    plot(t,xt(1,:)','-b','LineWidth', 2) 
    plot(t,xr(1,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Estado 1')
    box on; xlabel('t'); 
    grid on; box on;  hold off;
    
    subplot(4,1,2); grid on; hold on;
    plot(t,xp(2,:)','-r','LineWidth', 3) 
    plot(t,xt(2,:)','-b','LineWidth', 2) 
    plot(t,xr(2,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Estado 2');
    box on; xlabel('t');
    grid on; box on;  hold off;
    
    subplot(4,1,3); grid on; hold on;
    plot(t,xp(3,:)','-r','LineWidth', 3) 
    plot(t,xt(3,:)','-b','LineWidth', 2) 
    plot(t,xr(3,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Estado 3');
    grid on; box on;  hold off;
    
    subplot(4,1,4); grid on; hold on;
    plot(t,xp(4,:)','-r','LineWidth', 3) 
    plot(t,xt(4,:)','-b','LineWidth', 2) 
    plot(t,xr(4 ,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Estado 3');
    grid on; box on;  hold off;
    %---------------------------------Entrada------------------------------
    
    figure; subplot(2,1,1); grid on; hold on;
    plot(t,up(1,:)','-r','LineWidth', 3);
    plot(t,ut(1,:)','-b','LineWidth', 2);
    plot(t,ur(1,:)','--g','LineWidth', 1); 
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Entrada u_2');
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  hold off;
    
    subplot(2,1,2); grid on; hold on;
    plot(t,up(2,:)','-r','LineWidth', 3);
    plot(t,ut(2,:)','-b','LineWidth', 2);
    plot(t,ur(2,:)','--g','LineWidth', 1); 
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Entrada u_1');
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  hold off;

   
 