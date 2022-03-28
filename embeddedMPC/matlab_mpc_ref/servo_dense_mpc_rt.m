%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS - SERVOMOTOR
    clc; clear all; close all;
    
    % System information 
    K=39.08/27.92; 
    tau=1/27.92;
    T = 0.001;
    A = [exp(-T/tau),0; tau*(1-exp(-T/tau)),1];
    B = [K*(1-exp(-T/tau));T*K+tau*K*(exp(-T/tau)-1)];
    C = [0,1];
    x0= [0;0];
    
    % Constraints information
    N=4;                              Ns=10000;
    umin=-10;                         umax=10;
    xmin=[-5;-2];                     xmax=[5;2];
    Gamma=0.1; Omega=C'*C;  
    [Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);
    
    % Reference
    yref=[ones(1,Ns/4),-ones(1,Ns/4),ones(1,Ns/4),-ones(1,Ns/4)];
    
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
    G.linsol.solver='matlab';
    G.linsol.minres_iter=0;
    G.linsol.tol=1e-9;
    
    % QP solver information
    G.qpsol.solver='pdip';
    G.qpsol.pdip=20;
    
    % Save qps information
    G.sys.save='no';

    % MPC control
    [ut,xt,~]=general_dense_mpc_instate_rt(G);
    
    G.qpsol.solver='quadprog';
    [uj,xj,~]=general_dense_mpc_instate_rt(G);

    % LQR control with reference tracking
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

    % Plot input, state, and control signal
    t=1:1:Ns;
    % 
    figure; grid on; hold on;
    plot(t,xt(1,:)','r','LineWidth', 3) 
    plot(t,xj(1,:)','-b','LineWidth', 2) 
    plot(t,xr(1,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
    title('Estado 1 - Velocidad angular');
    box on; xlabel('t'); ylabel('wt'); 
    grid on; box on;  hold off;
    
    figure; grid on; hold on;
    plot(t,xt(2,:)','r','LineWidth', 3) 
    plot(t,xj(2,:)','-b','LineWidth', 2) 
    plot(t,xr(2,:),'--g','LineWidth', 1);
    legend('MPC_{pdip}','MPC_{quadprog}','LQR'); 
    title('Estado 2 - Posición angular')
    ylim([-2,2]);
    box on; xlabel('t'); ylabel('theta'); 
    grid on; box on;  hold off;
    %
    figure; grid on; hold on;
    plot(t,ut,'r','LineWidth', 3);
    plot(t,uj,'-b','LineWidth', 2);
    plot(t,ur,'--g','LineWidth', 1); 
    legend('MPC_{pdip}','MPC_{quadprog}','LQR'); 
    title('Entrada u_t (Vm)');
    ylim([-10,10]);
    ylabel('u_t'); xlabel('t'); 
    grid on; box on;  hold off;
    
