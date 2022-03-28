%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA
    % Dc-Motor example in: S. Lucia, D. Navarro, O. Lucia, P. Zometa, and R. Findeisen, 
    % “Optimized FPGA Implementation of Model Predictive Control for Embedded Systems Using 
    % High-Level Synthesis Tool,” IEEE Trans. Ind. Informatics, vol. 14,
    % no. 1, pp. 137–145, 2018.
    
    clc; clear all; close all;
    
    % Systems information
    ts=0.004; T=0.06; K=0.15;
    A=[1,ts;0,1-ts/T];
    B=[0;ts*(K/T)];
    x0= [0.8;-0.4];

    % Constraints information
    N=40;                             
    Ns=100;
    Omega=[1.1e4,0;0,2.9e1];          
    umin=-100;          
    umax=100;
    Gamma=2.4e-1;
    OmegaN = Omega;
    [Linf,~,~]=dlqr(A,B,Omega,Gamma);
    
    
    % Systems information
    G.sys.A=A;
    G.sys.B=B;
    G.sys.x0=x0;
    G.sys.Ns=Ns;

    
    % Cost function information
    G.cost.N=N;
    G.cost.OmegaN=OmegaN;
    G.cost.Omega=Omega;
    G.cost.Gamma=Gamma;
    
    % Constraints information
    G.cons.umin=umin;
    G.cons.umax=umax;
    
    % Linear solver information
    G.linsol.solver='lschol';
    G.linsol.minres_iter=30;
    G.linsol.tol=1e-9;
    
    % QP solver information
    G.qpsol.solver='pdip';
    G.qpsol.pdip=20;
    
    % Cálculo de la señal de control óptima y estado - control MPC
    [ut,xt,~]=general_dense_mpc_instate_rt(G);
    
   
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
