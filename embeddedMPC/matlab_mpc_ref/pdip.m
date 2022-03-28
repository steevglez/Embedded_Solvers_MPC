function [tk,val,tkt,iter]=pdip(H,h,Mx,cx,IT,tol,solver,varargin)
%===================================================================
% Quadratic programing (QP) problem solution:
%                   min 0.5 tk'*H*tk + h'*tk
%                   s.t     F*tk  = f
%                           G*tk <= g
% Arrange the equality and inequality constraints into a inequality 
%
% [tk,val,tkt,iter]=pdip(H,h,Mx,cx,IT,tol,mrmax,solver)
%                   min 0.5 tk'*H*tk + h'*tk
%                   s.t     Mx*tk <= cx 
%                           Mx=[G;F;-F] 
%                           cx=[g;f;-f]                       
%   tk     : Optimal solution of qp problem
%   val    : Value of the cost function in tk optimal
%   tkt    : Evolution of tk
%   iter   : Iterations of cgrad algorithm
%   IT     : Iterations of pdip algorithm
%   solver : Linear solver for Az=b
%   tol    : Tolerance of solution Az=b
%   mrmax  : Iterations of minres solver
%===================================================================

% -------- Primal-dual interior point algorithm --------------------
    [n,~]=size(H); 
    [i,~]=size(Mx);
%------------ Initial parameters -----------------------------------
    tk=1*ones(n,1);  lk=0.5*ones(i,1); sk=0.5*ones(i,1);     
    em=ones(i,1); iter=0; tkt(:,1)=tk; 
%------------- Iterations of algorithm -----------------------------
    sgk=0.5; zko=zeros(size(h,1),1);
	if nargin==8; mrmax=varargin{1}; end
	
    for k=1:1:IT
        % -------------- Build Ak ----------------------- 
        RK=diag(lk./sk);    
        RKI=diag(sk./lk);                    
        Ak=H+Mx'*RK*Mx;
        % -------------- Build bk -----------------------
        muk=(lk'*sk)/i;  % Den=T*l+2*p (MPC)              
        HK=-H*tk-h-Mx'*lk;                              
        GK=-Mx*tk+cx-sk; 
        TK=-lk.*sk+sgk*muk*em; 
        bk= HK+Mx'*RK*(GK-TK./lk);                             
        % ---------- Solve teh system Ak*zk=bk ----------
        if strcmp('minres',solver)
            [zk,itr]=minres(Ak,bk,zko,mrmax,tol);
            iter(1,k)=itr;
        end
        if strcmp('cgrad',solver)
            [zk,itr]=cgrad(Ak,bk,zko,size(Ak,1),tol);
            iter(1,k)=itr;
        end
        if strcmp('lschol',solver)
            zk=lschol(Ak,bk);
        end
        if strcmp('matlab',solver)
            zk=Ak\bk;
        end
        % ----- Calculate Delta_lk and Delta_sk----------
        Dtk=zk(1:n,1);
        Dlk=-RK*(GK-TK./lk-Mx*Dtk);
        Dsk=TK./lk-RKI*Dlk;
        % Find max ak in (0,1];
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
