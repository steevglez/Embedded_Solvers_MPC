# Model Predictive Control en C++ (CPU). 

Estos códigos C++ implementan MPC y OSQP.

## Requerimientos

Para compilar los codigos de OSQP se necesita tener instalado el programa **cmake**.

### Instalar cmake

- Se puede instalar **cmake** por medio del comando:

  `sudo apt-get -y instalar cmake`
  
### Script

Se incluye un script para la compilacion de OSQP y codigos que implementan MPC, ademas de la ejecucion de un ejemplo con N=4. Para ejecutar este script se deben otorgar permisos de ejecucion por medio de:

  `chmod +x setup_osqp.sh`
  
Luego, ejecutar el script por medio de:

`./setup_osqp.sh`

Esto tomara los archivos del solver, los compilara, eliminara algunos archivos para evitar conflictos entre los makefile, finalmente ejecutara el archivo makefile del codigo MPC y ejecutara a estos codigos a modo de ejemplo. La compilacion de archivos del solver y ejecucion de los codigos MPC se puede visualizar por medio de la terminal.

## Modo de uso - Ejecucion, muestras(samples) y salida de datos

Una vez ejecutado el script y verificar la compilacion y ejecucion del ejemplo por terminal, se obtendra el ejecutable **runner** con horizonte de prediccion N=4 por defecto. Para ejecutar este programa con distintas configuraciones de muestras, salidas de datos e iteraciones, se debe ejecutar con el siguiente comando

`./runner samplesMPC_N[HOR_PREDIC_MUESTRAS].txt [OUTPUT_DATOS] [OUTPUT_TIEMPO] [NUM_ITERACIONES] [NUMERO_EJECUCION]`

Donde:

samplesMPC_N[HOR_PREDIC_MUESTRAS].txt - 

[OUTPUT_DATOS] -

[OUTPUT_TIEMPO] -

[NUM_ITERACIONES] -

[NUMERO_EJECUCION] -

========================================

### Muestras(samples) y salida de datos

Cada ejecucion de codigo debe ir junto a una archivo de muestras, samplesMPC_N4.txt en este caso.

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
