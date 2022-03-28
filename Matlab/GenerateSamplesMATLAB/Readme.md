# Model Predictive Control en MATLAB. 

Estos códigos MATLAB se utilizan para tener una referencia del correcto funcionamiento de MPC, PDIP y 3 linear solvers. Esta referencia consiste en entradas y salidas esperadas para cada algoritmo, las se guardan en archivos .txt para que los códigos C++ usen las mismas entradas y comparen los resultados obtenidos.

## Requerimientos


Para utilizar estos códigos se necesita tener instalados estos dos toolbox: **Control System Toolbox** y **Optimization Toolbox**


## Modo de uso


Antes de correr el código *servoMPCReferenceTracking.m* se deben modificar los siguientes parámetros según se necesite: 
  - `N = [2, 3, 4, 5, 8, 10];` Los horizontes de predicción que se usar para MPC.
  - `saveMat.MPC = true;` Si se guardan las entradas y salidas para MPC en *samples/samplesMPC_N[N].txt*. 
  - `saveMat.PDIP = true;` Si se guardan las entradas y salidas para PDIP en *samples/samplesPDIP_N[N].txt*.
  - `saveMat.LS = true;` Si se guardan las entradas y salidas para LS en *samples/samplesLS_N[N].txt*. Esta opción en particular hace que el código sea ~5 veces más lento.
  - `compare = false;` Si se comparan los resultados de MPC con QuadProg y LQR.

Una vez se tenga la configuración lista se puede correr el código. 

## Descripción de cada código

+ *ABcal.m*
  - Calcula las matrices Acal y Bcal.
+ *cgrad.m*
  - Implementación 'Conjugate Gradient' para resolver sistemas de ecuaciones lineales.
+ *controlMPC.m*
  - Implementación de 'Model Predictive Control' para controlar un proceso.
+ *myChol.m*
  - Implementación 'Cholesky Decomposition' para resolver sistemas de ecuaciones lineales.
+ *myMinres.m*
  - Implementación 'Minimal Residual Method' para resolver sistemas de ecuaciones lineales.
+ *pdip.m*
  - Implementación de 'Primal Dual Interior Point' para resolver problemas de programación quadrática (QP).
+ *servoMPCReferenceTracking.m*
  - Código que simula el funcionamiento y control de un servo motor utilizando MPC.
+ *setup_mpc.m* 
  - Configuración inicial de MPC.
+ *stationaryStateValues* 
   - Calcula los valores en estado estacionario para las variables *x* y *u*.
+ *writeLSSamples.m* 
  - Genera un archivo .txt con las entradas y salidas esperadas para linear solvers.
+ *writeMPCSamples.m* 
  - Genera un archivo .txt con las entradas y salidas esperadas para PDIP.
+ *writePDIPSamples.m*
  - Genera un archivo .txt con las entradas y salidas esperadas para MPC.
