function [u,saveMat]=controlMPC(x,yref,xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Omg,H,Mx,N,qpSol,linSol,iterPDIP,iterMINRES, saveMat)
    % quadprog options
    options =  optimset('Display','off');

    % Change of variables
    [x_star,u_star] = stationaryStateValues(A,B,C,yref);
    
    a=ones(N,1)*(umin-u_star);                    
    b=ones(N,1)*(umax-u_star);
    x_tilde=x-x_star;
    h=(2*x_tilde'*Acal'*Omg*Ocal)';
    Xmin=kron(ones(N,1),xmin-x_star); 
    Xmax=kron(ones(N,1),xmax-x_star); 
    c=[(Xmax-Acal*x_tilde);-(Xmin-Acal*x_tilde)];
    % Convertir restricciones caja en Desigualdad
    cx=[c;b;-a];
    % Solve QP problem
    if strcmp(qpSol,'pdip')
        [u_tilde,~,~,~,saveMat]=pdip(H,h,Mx,cx,iterPDIP,iterMINRES,linSol,1e-9,saveMat);
    end
    if strcmp(qpSol,'quadprog')
        u_tilde=quadprog(H,h,Mx,cx,[],[],[],[],[],options);
    end
    
    % Change of variables
    u=u_tilde(1)+u_star;
    % ---------------------------------------------------------------------
    if saveMat.PDIP
        i = saveMat.i;
        saveMat.h(:,:,(i-1)+1) = [h];
        saveMat.cx(:,:,(i-1)+1) = [cx];
        saveMat.u_tilde(:,:,(i-1)+1) = [u_tilde];
    end
end



