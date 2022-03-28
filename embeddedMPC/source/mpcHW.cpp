#include"mpcHW.h"

void mpcHW(elem x[N_SYS], elem yref[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[M_SYS], elem umax[M_SYS], elem Acal[N_SYS*N_QP][N_SYS], elem AcalOmgOcal[N_SYS][N_QP], elem H[N_QP][N_QP], elem Mx[N_SYS*N_QP*2+2*N_QP][N_QP], elem L_last[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u[N_SYS]){

	elem cx[6*N_QP];
	elem h[N_QP];
//#pragma HLS ARRAY_PARTITION variable=cx complete
#pragma HLS ARRAY_PARTITION variable=x complete
#pragma HLS ARRAY_PARTITION variable=yref complete
#pragma HLS ARRAY_PARTITION variable=l_last complete
	elem x_star[N_SYS];
	elem u_star[M_SYS];
	// [x_star,u_star] = Infinity(A,B,C,yref);
	set_x_star:
	for (int i=0; i<N_SYS; i++){
#pragma HLS UNROLL
		x_star[i] = L_last[i]*yref[0];
	}
	set_u_star:
	for (int i=0; i<M_SYS; i++){
#pragma HLS UNROLL
		u_star[i] = L_last[i+N_SYS]*yref[0];
	}

	elem a[N_QP];
	elem b[N_QP];
#pragma HLS ARRAY_PARTITION variable=a complete
#pragma HLS ARRAY_PARTITION variable=b complete
	//a=ones(N,1)*(umin-uinfy);
	//b=ones(N,1)*(umax-uinfy);
	set_a:
	for (int i=0; i<N_QP; i++){
#pragma HLS UNROLL
		a[i] = umin[0] - u_star[0];
	}
	set_b:
	for (int i=0; i<N_QP; i++){
#pragma HLS UNROLL
		b[i] = umax[0] - u_star[0];
	}
	elem x_tilde[N_SYS];
	//xnau=x-xinfy;
	set_x_tilde:
	for (int i=0; i<N_SYS; i++){
#pragma HLS UNROLL
		x_tilde[i] = x[i] - x_star[i];
	}

	//h=(2*xnau'*Acal'*Omg*Ocal)';
	get_h(AcalOmgOcal, x_tilde, h);

	//elem cx[6*N_QP];
	// c = [(Xmax-Acal*xnau);-(Xmin-Acal*xnau)]
	// cx = [c;b;a]
	get_cx(Acal, xmax, xmin, x_star, x_tilde, a, b, cx);
	elem tol = 1e-9;
	elem u_tilde[N_QP];
	pdipHW(H, h, Mx, cx, u_tilde, iterPDIP, iterLS, tol);
	u[0]=u_tilde[0]+u_star[0];
}

void get_h(elem AcalQOcal[N_SYS][N_QP], elem x_tilde[N_SYS], elem h[N_QP]){
#pragma HLS INLINE off
	set_h:
	for (int i=0; i<N_QP; i++){
#pragma HLS pipeline
		h[i] = 0;
		for (int j=0; j<N_SYS; j++){
			h[i] += 2*x_tilde[j]*AcalQOcal[j][i];
		}
	}
}
void get_cx(elem Acal[N_SYS*N_QP][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_star[N_SYS], elem x_tilde[N_SYS], elem a[N_QP], elem b[N_QP], elem cx[6*N_QP]){
	set_cx_c1:
	for (int i=0; i<2*N_QP; i++){
#pragma HLS pipeline
		cx[i] =   xmax[i%N_SYS]- x_star[i%N_SYS];
		for (int c=0; c<N_SYS; c++){
			cx[i] -= Acal[i%(N_SYS*N_QP)][c]*x_tilde[c];
		}
	}
	set_cx_c2:
	for (int i=0; i<2*N_QP; i++){
#pragma HLS pipeline
		cx[i+2*N_QP] =   - (xmin[i%N_SYS]- x_star[i%N_SYS]);
		for (int c=0; c<N_SYS; c++){
			cx[i+2*N_QP] += Acal[i%(N_SYS*N_QP)][c]*x_tilde[c];
		}
	}
	set_cx_b:
	for (int i=0; i<N_QP; i++){
#pragma HLS pipeline
		cx[i+4*N_QP] = b[i];
	}
	set_cx_a:
	for (int i=0; i<N_QP; i++){
#pragma HLS pipeline
		cx[i+5*N_QP] = -a[i];
	}
}
