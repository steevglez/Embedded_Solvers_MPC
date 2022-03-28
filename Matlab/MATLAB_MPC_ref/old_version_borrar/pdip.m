function [tk,val,tkt,iter]=pdip(H,h,Mx,cx,IT,tol,mrmax,solver)
%==========================================================================
% Solución del problema de programción cuadrática:
%                   min 0.5 tk'*H*tk + h'*tk
%                   s.t     F*tk  = f
%                           G*tk <= g
% Para la solución se agrupan todas las restricciones en una desigualdad
%
% [tk,val,tkt,iter]=pdip(H,h,Mx,cx,IT,tol,mrmax,solver)
%                   min 0.5 tk'*H*tk + h'*tk
%                   s.t     Mx*tk <= cx con 
%                           Mx=[G;F;-F] 
%                           cx=[g;f;-f]                       
%   tk     : Solución optima del problema de programción cuadrática
%   val    : Valor del funcional en la solución optima tk
%   tkt    : Evolución de la solución hasta llegar a ser optima
%   iter   : Número de iteraciones internas del algoritmo cg_wip
%   IT     : Número máximo de iteraciones del algoritmo de punto interior
%   solver : Solver para el sistema de equaciones Az=b
%   tol    : Tolerancia aceptada para la solución del sistema Az=b
%   mrmax  : Si el solver es minres, mrmax indica el máximo de iteraciones
%          : si el solver no es minres, el dato mrmax no se toma en cuenta
%          : sin embargo debe ser pasado a la función.
%==========================================================================

% -------- Algoritmo de Punto Interior --------------------------------
    [n,~]=size(H); 
    [i,~]=size(Mx);
%-------------Parámetros inicales--------------------------------------
    tk=1*ones(n,1);  lk=0.5*ones(i,1); sk=0.5*ones(i,1);     
    em=ones(i,1); iter=0; tkt(:,1)=tk; 
%-------------Iteraciones del algoritmo--------------------------------
    sgk=0.5; zko=zeros(size(h,1),1);
    for k=1:1:IT
        % ------------Construir Ak----------------------
        %Lk=diag(lk);                             
        %Sk=diag(sk); 
        %RK=Lk*Sk^-1; 
        RK=diag(lk./sk);    
        RKI=diag(sk./lk);                    
        Ak=H+Mx'*RK*Mx
        % ------------Construir bk----------------------
        muk=(lk'*sk)/i;  % Den=T*l+2*p (MPC)              
        HK=-H*tk-h-Mx'*lk;                              
        GK=-Mx*tk+cx-sk; 
        TK=-lk.*sk+sgk*muk*em; 
        bk= HK+Mx'*RK*(GK-TK./lk);                             
        % ----------Resolver el sistema Ak*zk=bk---------
        if strcmp('minres',solver)
            [zk,itr]=min_res(Ak,bk,zko,mrmax,tol);
            iter(1,k)=itr;
        end
        if strcmp('cg_wip',solver)
            [zk,itr]=cg_wip(Ak,bk,zko,size(Ak,1),tol);
            iter(1,k)=itr;
        end
        if strcmp('lschol',solver)
            zk=lschol(Ak,bk);
        end
        if strcmp('matlab',solver)
            zk=Ak\bk;
        end
        % ---------Cacular Delta_lk y  Delta_sk----------
        Dtk=zk(1:n,1);
        Dlk=-RK*(GK-TK./lk-Mx*Dtk);
        Dsk=TK./lk-RKI*Dlk;
        % Encontrar max ak en (0,1];
        bt=0.99999;
        alp_lk=1;
        alp_sk=1;
        idx=find(Dlk<0);
        if (isempty(idx)==0)
            alp_lk = bt/(max(alp_lk,max(-Dlk(idx)./lk(idx))));
        end
        idx=find(Dsk<0);
        if (isempty(idx)==0)
            alp_sk = bt/(max(alp_sk,max(-Dsk(idx)./sk(idx))));
        end
        alp=min(alp_lk,alp_sk);
        tk=tk+alp*Dtk;    
        lk=lk+alp*Dlk;   
        sk=sk+alp*Dsk;
        tkt(:,k+1)=tk; zko=zk;
    end
    val=0.5*tk'*H*tk+h'*tk;
end   
