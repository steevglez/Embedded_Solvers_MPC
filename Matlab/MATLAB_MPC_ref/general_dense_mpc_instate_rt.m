function [uopt,xt,yt]=general_dense_mpc_instate_rt(G)
    % Model Predictive control for input and state constraints and
    % reference tracking. General version. Angel Cedeño Nieto
    
    % Quadprog configuration
    options =  optimset('Display','off');
    
    % Systems information
    A=G.sys.A;
    B=G.sys.B;
    x0=G.sys.x0;
    Ns=G.sys.Ns;
    
    n=size(A,1);
    m=size(B,2);
    
    if isfield(G.sys,'C')
        C=G.sys.C;
    else
        C=zeros(1,n);  
    end
    p=size(C,1);
    
    % Cost function information
    N=G.cost.N;
    OmegaN=G.cost.OmegaN;
    Omega=G.cost.Omega;
    Gamma=G.cost.Gamma;

    
    % Linear solver information
    linsol=G.linsol.solver;
    tol=G.linsol.tol;
    
    % QP solver information
    qpsol=G.qpsol.solver;
    qpiter=G.qpsol.pdip;
    
    % Extended system matrices Acal, Bcal for t=1...N
    [Acal,AN,Bcal,BN] = cal_matrices_extended_system(A,B,N);
    
    % Cost function matrices - qp problem
    Gmm=kron(eye(N),Gamma);
    Omg=blkdiag(kron(eye(N-1),Omega),OmegaN);
    
    % H: matrix of cost Vn(x_0,vec{u})
    H=2*Gmm+2*Bcal'*Omg*Bcal;              
    H=(H+H')/2;
    
    % Defining the kind of constraints
    cons1=[0,0,0]; cons2=[0,0,0]; cons3=[0,0,0]; 
    if isfield(G.cons,'Nxmin') && isfield(G.cons,'Nxmax')
        Nxmin=G.cons.Nxmin;
        Nxmax=G.cons.Nxmax;
        cons1=[1,0,0]; 
    end
    if isfield(G.cons,'xmin') && isfield(G.cons,'xmax')
        xmin=G.cons.xmin;
        xmax=G.cons.xmax;
        cons2=[0,1,0];
    end
    if isfield(G.cons,'umin') && isfield(G.cons,'umax')
        umin=G.cons.umin;
        umax=G.cons.umax;
        cons3=[0,0,1]; 
    end
    cons=cons1+cons2+cons3;
    cass=cons*[4;2;1];
    
    % Mx: matrix of linear constraints Mx*u<=cx
    switch cass
        %only_finalstate_constraints
        case 4
            Mx=[BN;-BN]; 
        % only_state_constraints    
        case 2
            Mx=[Bcal;-Bcal];  
        % only_input_constraints    eye(N*m);-eye(N*m) 
        case 1
            Mx=[eye(N*m);-eye(N*m)];  
        % finalstate_state_constraints
        case 6
            Mx=[BN;-BN;Bcal;-Bcal];  
        % state_input_constraints
        case 3
            Mx=[Bcal;-Bcal;eye(N*m);-eye(N*m)];  
        % finalstate_input_constraints
        case 5
            Mx=[BN;-BN;eye(N*m);-eye(N*m)];
        % finalstate_state_input_constraints
        case 7
            Mx=[BN;-BN;Bcal;-Bcal;eye(N*m);-eye(N*m)]; 
    end
    
    % New x0 for reference tracking
    xt(:,1)=x0;
    
    % Control loop
    uopt=zeros(m,Ns);
    yt=zeros(p,Ns);
    
    for i=1:1:Ns
        
        % Reference tracking?
        if isfield(G.sys,'yref')
            yref=G.sys.yref;
            % Get the infinity value
            [xinfy,uinfy] = stationary_state_values_rt(A,B,C,zeros(m,1),yref(:,i));
        else
            xinfy=zeros(n,1); uinfy=zeros(m,1); 
        end
 
        % Initial conditions
        x0nau=x0-xinfy; 

        % f: vector of cost vn(x_0,vec{u})
        h=(2*x0nau'*Acal'*Omg*Bcal)';
        
        % cx: vector of linear constraints Mxu<=cx
        switch cass
            %only_finalstate_constraints
            case 4
                aN=((Nxmin-xinfy)-AN*x0nau); 
                bN=((Nxmax-xinfy)-AN*x0nau);
                cx=[bN;-aN];
            % only_state_constraints  
            case 2
                Xmin=kron(ones(N,1),xmin-xinfy); 
                Xmax=kron(ones(N,1),xmax-xinfy);
                ax=(Xmin-Acal*x0nau);
                bx=(Xmax-Acal*x0nau);
                cx=[bx;-ax];
            % only_input_constraints  
            case 1
                a=kron(ones(N,1),umin-uinfy);                   
                b=kron(ones(N,1),umax-uinfy);
                cx=[b;-a];
            % finalstate_state_constraints
            case 6
                Xmin=kron(ones(N,1),xmin-xinfy); 
                Xmax=kron(ones(N,1),xmax-xinfy);
                aN=((Nxmin-xinfy)-AN*x0nau);
                bN=((Nxmax-xinfy)-AN*x0nau);
                ax=(Xmin-Acal*x0nau);
                bx=(Xmax-Acal*x0nau);
                cx=[bN;-aN;bx;-ax];
            % state_input_constraints
            case 3
                Xmin=kron(ones(N,1),xmin-xinfy); 
                Xmax=kron(ones(N,1),xmax-xinfy);
                ax=(Xmin-Acal*x0nau);
                bx=(Xmax-Acal*x0nau);
                a=kron(ones(N,1),umin-uinfy);                    
                b=kron(ones(N,1),umax-uinfy);
                cx=[bx;-ax;b;-a];
            % finalstate_input_constraints
            case 5
                aN=((Nxmin-xinfy)-AN*x0nau); 
                bN=((Nxmax-xinfy)-AN*x0nau);
                a=kron(ones(N,1),umin-uinfy);                    
                b=kron(ones(N,1),umax-uinfy);
                cx=[bN;-aN;b;-a];
            % finalstate_state_input_constraints
            case 7
                Xmin=kron(ones(N,1),xmin-xinfy); 
                Xmax=kron(ones(N,1),xmax-xinfy);
                aN=((Nxmin-xinfy)-AN*x0nau);
                bN=((Nxmax-xinfy)-AN*x0nau);
                ax=(Xmin-Acal*x0nau);
                bx=(Xmax-Acal*x0nau);
                a=kron(ones(N,1),umin-uinfy);                    
                b=kron(ones(N,1),umax-uinfy);
                cx=[bN;-aN;bx;-ax ;b;-a]; 
        end

        % Solving qp problem
        if strcmp(qpsol,'pdip')
            if strcmp(linsol,'minres')
                minres_iter=G.linsol.minres_iter;
                unau=pdip(H,h,Mx,cx,qpiter,tol,linsol,minres_iter);
            else
                unau=pdip(H,h,Mx,cx,qpiter,tol,linsol);
            end
        end
        if strcmp(qpsol,'quadprog')
            unau=quadprog(H,h,Mx,cx,[],[],[],[],[],options);
        end
        ut=unau(1:m)+uinfy;
        uopt(:,i)=ut;
        
        % Initial state for the next iteration
        xt(:,i+1)=A*xt(:,i)+B*ut;
        yt(:,i)=C*xt(:,i);
        x0=xt(:,i+1); 
        
        % Saving qp problems
        if isfield(G.sys,'save')
            if strcmp(G.sys.save,'yes')
                ht(:,:,i)=h;
                cxt(:,:,i)=cx;
                sqp(:,:,i)=unau;
            end
        end
    end
    
    % Saving qp problems
    if isfield(G.sys,'save')
        if strcmp(G.sys.save,'yes')
            save(['H_',int2str(N)],'H');
            save(['ht_',int2str(N)],'ht');
            save(['Mx_',int2str(N)],'Mx');
            save(['cx_',int2str(N)],'cxt');
            save(['sqp_',int2str(N)],'sqp');
        end
    end

    xt=xt(:,1:end-1);
end

