function [tk,val,tkt,iter,saveMat]=pdip(H,h,M,c,iterPDIP,iterMINRES,solver,tol,saveMat)
%===================================================================
% [tk,val,tkt,iter]=pdip(H,h,M,c,IT,tol,solver, mrmax)
%                   min 0.5 tk'*H*tk + h'*tk
%                   s.t     M*tk <= c                       
%   tk      : Optimal solution of qp problem
%   val     : Value of the cost function in tk optimal
%   tkt     : Evolution of tk
%   iterPDIP : Iterations of pdip algorithm
%   solver  : Linear solver for Az=b
%   tol     : Tolerance of solution Az=b
%   iterLS  : Max iterations for minres linear solver
%   saveMat : Options to save Mat (saveMat.LS saveMat.outputMat)
%===================================================================

% -------- Primal-dual interior point algorithm --------------------
    [n,~]=size(H); 
    [i,~]=size(M);
%------------ Initial parameters -----------------------------------
    tk=1*ones(n,1);
    lk=0.5*ones(i,1);
    sk=0.5*ones(i,1);
    em=ones(i,1);
    tkt(:,1)=tk; 
%------------- Iterations of algorithm -----------------------------
    sgk=0.5;
    zko=zeros(n,1);
    
%------------- Outputs for reference  ------------------------------
if saveMat.LS
    Ak_out = [];
    bk_out = [];
    zk_out = [];
end

    for iter=1:1:iterPDIP
        % -------------- Build Ak ----------------------- 
        RK=diag(lk./sk);    
        RKI=diag(sk./lk);                   
        Ak=H+M'*RK*M;
        % -------------- Build bk -----------------------
        muk=(lk'*sk)/i;  % Den=T*l+2*p (MPC)              
        HK=-H*tk-h-M'*lk;                              
        GK=-M*tk+c-sk; 
        TK=-lk.*sk+sgk*muk*em;
        bk= HK+M'*RK*(GK-TK./lk);  
        % ---------- Solve the system Ak*zk=bk ----------
        if strcmp('minres',solver)
            [zk,~]=myMinres(Ak,bk,zko,iterMINRES,tol);
        end
        if strcmp('cgrad',solver)
            [zk,~]=cgrad(Ak,bk,zko,size(Ak,1),tol);
        end
        if strcmp('chol',solver)
            zk=myChol(Ak,bk);
        end
        if strcmp('matlab',solver)
            zk=Ak\bk;
        end
        % ----- Calculate Delta_lk and Delta_sk----------
        Dlk=-RK*(GK-TK./lk-M*zk);
        Dsk=TK./lk-RKI*Dlk;
        % ----------- Find best alp ---------------------
        bt=0.99999;
        alp_lk=1;
        alp_sk=1;
        idx=find(Dlk<0);
        if (isempty(idx)==0)
            alp_lk = max(alp_lk,max(-Dlk(idx)./lk(idx)));
        end
        idx=find(Dsk<0);
        if (isempty(idx)==0)
            alp_sk = max(alp_sk,max(-Dsk(idx)./sk(idx)));
        end
        alp=bt/max(alp_lk,alp_sk);
        
        % ----------- Prepare the next iteration ---------
        tk=tk+alp*zk;    
        lk=lk+alp*Dlk; 
        sk=sk+alp*Dsk;
        zko=zk;
        tkt(:,iter+1)=tk; 
        % update reference outputs
        if saveMat.LS
            Ak_out = cat(3, Ak, Ak_out);
            bk_out = cat(3, bk, bk_out);
            zk_out = cat(3, zk, zk_out);
        end
    end
    val=0.5*tk'*H*tk+h'*tk;
    
    if saveMat.LS
        saveMat.Ak = cat(3, saveMat.Ak, Ak_out);
        saveMat.bk = cat(3, saveMat.bk, bk_out);
        saveMat.zk = cat(3, saveMat.zk, zk_out);
    end
end   