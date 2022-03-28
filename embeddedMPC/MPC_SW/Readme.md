# Model Predictive Control en C++ (CPU). 

Estos códigos C++ implementan MPC y PDIP. PDIP hace uso de OpenBLAS para algunas operaciones de matrices y dos linear solvers.

## Requerimientos


Para utilizar estos códigos se necesita tener la biblioteca **OpenBLAS** instalada.

### Instalar OpenBLAS

- Clonar el repositorio de OpenBLAS

  `git clone git@github.com:xianyi/OpenBLAS.git`
- Compilar OpenBLAS

  `cd OpenBLAS`

  `make`

  `make install`

  por defecto queda instalado en */opt/OpenBLAS*

### Muestras (samples)
Se tienen que tener las muestras generadas por el código *servoMPCReferenceTracking.m* con el horizonte de predicción que se va a utilizar. Para generar estas muestras seguir las instrucciones en GenerateSamplesMATLAB/Readme.md

## Modo de uso

Antes de compilarse debe modificar en el header file *source/specs.h* el horizonte de predicción que se utilizará `#define N_QP 4` y se debe des-comentar la linea `#define ELEM_FLOAT` o `#define ELEM_DOUBLE` para fijar el tamaño de punto flotante que se va a ocupar.

Para compilar solo se debe correr el makefile:

 `make`
 
 Si la biblioteca OpenBLAS no fue instalada en */opt/OpenBLAS* se debe modificar las lineas de *LIBS* y *INCLUDES* del makefile.

**Eventualmente el makefile va a tener la opción de compilar para ARM**

Para correr el código se debe indicar el .txt con las muestras, por ejemplo para un horizonte N=4:

`./runner ../GenerateSamplesMATLAB/samples/samplesMPC_N4.txt`


## Descripción de cada código

TO DO
<!---
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
  - Genera un archivo .txt con las entradas y salidas esperadas para MPC.}

-->
