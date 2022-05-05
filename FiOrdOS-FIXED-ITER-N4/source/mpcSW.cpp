#include"mpcSW.h"
#include <string.h>

void mpcSW(elem x[N_SYS], elem yref[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[M_SYS], elem umax[M_SYS], elem A[N_SYS][N_SYS], elem B[N_SYS], elem Acal[N_SYS*N_QP][N_SYS], elem AcalOmgOcal[N_SYS][N_QP], elem H[N_QP *N_QP], elem Mx[N_SYS*N_QP*2+2*N_QP *N_QP], elem L_last[N_SYS+P_SYS], elem u[N_SYS], fiordos4x4_Cparam_Params *params, fiordos4x4_Cparam_Settings *settings, fiordos4x4_Cparam_Result *result, fiordos4x4_Cparam_Work *work, int *started, double output_array[], tictoc *SOLVERTimer, int sample){

	elem cx[6*N_QP];
	elem h[N_QP];
	elem x_star[N_SYS];
	elem u_star[M_SYS];
	elem a[N_QP];
	elem b[N_QP];
	elem x_tilde[N_SYS];
	elem u_tilde[N_QP];
	
	// Calculo de seguimiento de referencia distinta de cero
	// MATLAB: [x_star,u_star] = Infinity(A,B,C,yref);
	for (int i=0; i<N_SYS; i++){
		x_star[i] = L_last[i]*yref[0];
	}
	for (int i=0; i<M_SYS; i++){
		u_star[i] = L_last[i+N_SYS]*yref[0];
	}
	// Calculo de restriccion en estado estacionario a
	// MATLAB: a=ones(N,1)*(umin-uinfy);
	for (int i=0; i<N_QP; i++){
		a[i] = umin[0] - u_star[0];
	}
	// Calculo de restriccion en estado estacionario b
	// MATLAB: b=ones(N,1)*(umax-uinfy);	
	for (int i=0; i<N_QP; i++){
		b[i] = umax[0] - u_star[0];
	}
	// Cambio de variable para estado estacionario 
	// MATLAB: x_tilde=x-x_star;
	for (int i=0; i<N_SYS; i++){
		x_tilde[i] = x[i] - x_star[i];
	}
	
	// Obtener vector funcional o lineal h
	// MATLAB: h=(2*x_tilde'*Acal'*Omg*Ocal)';
	// AclaOmgOcal = Acal'*Omg*Ocal
    	get_h(AcalOmgOcal, x_tilde, h);

	// Obtener vector de restricciones c
	// MATLAB: c = [(Xmax-Acal*x_tilde);-(Xmin-Acal*x_tilde)];
	// MATLAB: cx = [c;b;a];
	get_cx(Acal, xmax, xmin, x_star, x_tilde, a, b, cx);

	// Llamada a funcion del que toma las estructuras de datos y las pasa al solver
	// MATLAB: [u_tilde,~,~,~,saveMat]=pdip(H,h,Mx,cx,iterPDIP,1e-9,linSol,mrmax,saveMat); - PDIP REEMPLAZADO POR EL SOLVER FIORDOS
	fiordosSW(H, h, Mx, cx, u_tilde, params, settings, result, work, started, SOLVERTimer);

	// Desnormalizacion de actuacion
	// MATLAB: u=u_tilde(1)+u_star;
	u[0]=u_tilde[0]+u_star[0];
	
	// Almacenamiento de datos
	output_array[3*sample] = u[0];
	output_array[3*sample+1] = x[0];
	output_array[3*sample+2] = x[1];
	
	// Actualizacion de estados dada la actuacion calculada
	x[0] = A[0][0]*x[0] + A[0][1]*x[1] + B[0]*u[0];
	x[1] = A[1][0]*x[0] + A[1][1]*x[1] + B[1]*u[0];
}

void get_h(elem AcalQOcal[N_SYS][N_QP], elem x_tilde[N_SYS], elem h[N_QP]){
	for (int i=0; i<N_QP; i++){
        h[i] = 0;
		for (int j=0; j<N_SYS; j++){
            h[i] += 2*x_tilde[j]*AcalQOcal[j][i];
		}
	}
}
void get_cx(elem Acal[N_SYS*N_QP][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_star[N_SYS], elem x_tilde[N_SYS], elem a[N_QP], elem b[N_QP], elem cx[6*N_QP]){
	for (int i=0; i<2*N_QP; i++){
		cx[i] =   xmax[i%N_SYS]- x_star[i%N_SYS];
		for (int c=0; c<N_SYS; c++){
			cx[i] -= Acal[i%(N_SYS*N_QP)][c]*x_tilde[c];
		}
	}
	
	for (int i=0; i<2*N_QP; i++){
		cx[i+2*N_QP] =   - (xmin[i%N_SYS]- x_star[i%N_SYS]);
		for (int c=0; c<N_SYS; c++){
			cx[i+2*N_QP] += Acal[i%(N_SYS*N_QP)][c]*x_tilde[c];
		}
	}
	
	for (int i=0; i<N_QP; i++){
		cx[i+4*N_QP] = b[i];
	}
	
	for (int i=0; i<N_QP; i++){
		cx[i+5*N_QP] = -a[i];
	}
}
