%% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
    % Chain of masses example in: S. Lucia, D. Navarro, O. Lucia, P. Zometa, and R. Findeisen, 
    % “Optimized FPGA Implementation of Model Predictive Control for Embedded Systems Using 
    % High-Level Synthesis Tool,” IEEE Trans. Ind. Informatics, vol. 14,
    % no. 1, pp. 137–145, 2018.


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
    N=10;                             Ns=150;                                  
    umin=-0.5*ones(3,1);              umax=0.5*ones(3,1);
    xmin=-4*ones(12,1);               xmax=4*ones(12,1);
    Nxmin=xmin;                       Nxmax=xmax;
    Gamma=eye(3);                     Omega=eye(12);           OmegaN=Omega;
    
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
    G.cons.xmin=xmin;
    G.cons.xmax=xmax;
    G.cons.umin=umin;
    G.cons.umax=umax;
    G.cons.Nxmin=Nxmin;
    G.cons.Nxmax=Nxmax;
    
    % Linear solver information
    G.linsol.solver='cgrad';
    G.linsol.minres_iter=30;
    G.linsol.tol=1e-9;
    
    % QP solver information
    G.qpsol.solver='pdip';
    G.qpsol.pdip=20;
       
    % Cálculo de la señal de control óptima y estado - control MPC
    [ut,xt,~]=general_dense_mpc_instate_rt(G);
    
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