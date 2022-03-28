#include "pdipSW.h"
#include "utils.h"

void pdipSW  (elem H[N_QP*N_QP], elem h[N_QP], elem M[I_QP*N_QP], elem c[I_QP], elem tk[N_QP], int iterQp, int iterLs, elem tol, int solver){
    
    elem sgk = 0.5;
    elem bt = 0.99999;
    elem alp = 1.0;

    elem bk[N_QP];
    elem muk;
    elem HK[N_QP];
    elem HK_aux1[N_QP];
    elem HK_aux2[N_QP];
    elem GK[I_QP];
    elem GK_aux[I_QP];
    elem RK[I_QP*I_QP] = {0};
    elem Ak[N_QP*N_QP];
    elem Ak_aux[N_QP*I_QP];
    elem bk_aux1[I_QP];
    elem zk[N_QP];
    elem zko[N_QP] = {0};
    elem Dlk[I_QP];
    elem Dlk_aux1[I_QP];
    elem Dlk_aux2[I_QP];
    elem Dsk[I_QP];

    elem lk[I_QP];
    elem sk[I_QP];

    // MATLAB: tk=1*ones(n,1);
    for (int i=0; i<N_QP; i++){
        tk[i] = 1;
    }
    // MATLAB: lk=0.5*ones(i,1); 
    // MATLAB: sk=0.5*ones(i,1);
    for (int i=0; i<I_QP; i++){
        lk[i] = .5;
        sk[i] = .5;
    }

    // PDIP iterations
    for (int iter=0; iter<iterQp; iter++){
        // -------------- Build Ak -----------------------
        
        // MATLAB: RK=diag(lk./sk);
        // RK is already inilialized with ceros, then its only 
        // necesary to fill the diagonal
        for (int r=0; r<I_QP; r++){
            RK[r*I_QP+r] = lk[r]/sk[r];

        }
        // MATLAB: Ak=H+M'*RK*M; 
        // Ak_aux = M'*RK
        mmultSW(M, CblasTrans, RK, CblasNoTrans, Ak_aux, N_QP, I_QP, I_QP);
        // Ak = Ak_aux*M
        mmultSW(Ak_aux, CblasNoTrans, M, CblasNoTrans, Ak, N_QP, I_QP, N_QP);
        // Ak = Ak+H
        AXPY(N_QP*N_QP, 1, H, 1, Ak, 1);

        // -------------- Build bk -----------------------

        // MATLAB: muk=(lk'*sk)/i;
        // muk = (lk'*sk)
        muk = DOT(I_QP, lk, 1, sk, 1);
        // muk = muk/i
        muk = muk/I_QP;

        // MATLAB: HK=-H*tk-h-M'*lk;
        // HK_aux1 = H*tk
        mmultSW(H, CblasNoTrans, tk, CblasNoTrans, HK_aux1, N_QP, N_QP, 1);
        // HK_aux2 = M'*lk
        mmultSW(M, CblasTrans, lk, CblasNoTrans, HK_aux2, N_QP, I_QP, 1);
        // HK = -HK_aux1 -h -HK_aux2
        for (int i=0; i<N_QP; i++){
            HK[i] = -HK_aux1[i] - h[i] - HK_aux2[i];
        }

        // MATLAB: GK=-M*tk+c-sk;
        // GK_aux = M*tk
        mmultSW(M, CblasNoTrans, tk, CblasNoTrans, GK_aux, I_QP, N_QP, 1);
        // GK = -GK_aux +c -sk
        for (int i=0; i<I_QP; i++){
            GK[i] = -GK_aux[i] + c[i] - sk[i];
        }

        // MATLAB: TK=-lk.*sk+sgk*muk*em;  em is an identity matrix
        // TK = -lk.*sk + sgk*muk
        elem TK[I_QP];
        for (int i=0; i<I_QP; i++){
            TK[i] = -lk[i]*sk[i] + sgk*muk;
        }

        // MATLAB: bk= HK+M'*RK*(GK-TK./lk); 
        // bk_aux1 = (GK-TK./lk)
        for (int i=0; i<I_QP; i++){
            bk_aux1[i] = GK[i]-(TK[i]/lk[i]);
        }
        // Ak_aux = M'*RK
        // bk = bk_aux2*bk_aux1
        mmultSW(Ak_aux, CblasNoTrans, bk_aux1, CblasNoTrans, bk, N_QP, I_QP, 1);
        // bk = Hk+bk
        AXPY(N_QP, 1, HK, 1, bk, 1);

        // --------- Solve the system Ak*zk=bk ----------

        if (solver==MINRES_SW){
            minresSW(Ak, bk, zko, zk, iterLs, tol);
        }
        else if (solver==CGRAD_SW){
            cgradSW(Ak, bk, zko, zk, tol);
        }
        else if (solver==CHOL_SW){
            cholSW(Ak, bk, zk);
        }
        else if (solver==LAPACKCHOL_SW){
            lapackcholSW(Ak, bk, zk);
        }
        else if (solver==LAPACKGESV_SW){
            lapackgesvSW(Ak, bk, zk);
        }
        else{
            cerr << "ERROR, wrong solver: " << solver << endl;
            exit(1);
        }

        // ----- Calculate Delta_lk and Delta_sk----------

        // MATLAB: Dlk=-RK*(GK-TK./lk-M*zk);
        // Dlk_aux1 = M*zk
        mmultSW(M, CblasNoTrans, zk, CblasNoTrans, Dlk_aux1, I_QP, N_QP, 1);
        // Dlk_aux2 = -GK + TK/lk + Dlk_aux1
        for (int i=0; i<I_QP; i++){
            Dlk_aux2[i] = (-GK[i] + TK[i]/lk[i] + Dlk_aux1[i]);
        }
        // Dlk = RK*Dlk_aux2
        mmultSW(RK, CblasNoTrans, Dlk_aux2, CblasNoTrans, Dlk, I_QP, I_QP, 1);

        // MATLAB: Dsk=TK./lk-RKI*Dlk;
        // Dsk = TK/lk - Dlk/RK_diag
        // RK_diag = RK(i,i) 
        for (int i=0; i<I_QP; i++){
            Dsk[i] =  TK[i]/lk[i] -Dlk[i]/RK[i*I_QP+i];
        }

        // ----------- Find best alp ---------------------

        alp = 1.0;
        for (int i=0; i<I_QP; i++){
            // MATLAB: idx=find(Dlk<0);
            // set Dlk_i to 0 if Dlk[i]>0. This way
            // a Dlk[i]>0 wont change the results. (same for Dsk_i)
            elem Dlk_i = (Dlk[i]<0) ? Dlk[i] : 0;
            // MATLAB: idx=find(Dsk<0);
            elem Dsk_i = (Dsk[i]<0) ? Dsk[i] : 0;

            // MATLAB: alp_lk = max(alp_lk,max(-Dlk(idx)./lk(idx)));
            // MATLAB: alp_sk = max(alp_sk,max(-Dsk(idx)./sk(idx)));
            // Searching for the best alp instead of searching for 
            // the best alp_sk and alp_lk separately. 
        	Dlk_i = -Dlk_i/lk[i];
        	Dsk_i = -Dsk_i/sk[i];

        	alp = (alp>Dlk_i)? alp : Dlk_i;
        	alp = (alp>Dsk_i)? alp : Dsk_i;
        }
        // MATLAB: alp=bt/max(alp_lk,alp_sk);
        // Here alp is already max(alp_lk,alp_sk)
        alp = bt/alp;

        // ----------- Prepare the next iteration ---------

        // MATLAB: tk=tk+alp*zk;
        AXPY(N_QP, alp, zk, 1, tk, 1);
        // MATLAB: zko=zk;
        COPY(N_QP, zk, zko);
        // MATLAB: lk=lk+alp*Dlk;
        AXPY(I_QP, alp, Dlk, 1, lk, 1);
        // MATLAB: sk=sk+alp*Dsk;
        AXPY(I_QP, alp, Dsk, 1, sk, 1);
    }
}

