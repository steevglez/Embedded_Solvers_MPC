# Model Predictive Control en C++ (CPU). 

Estos códigos C++ implementan MPC y OSQP.

## Requerimientos

Para compilar los codigos de OSQP se necesita tener instalado el programa **cmake**.

### Instalar cmake

- Se instala **cmake** por medio del comando:

  `sudo apt-get -y install cmake`
  
### Uso de script de setup

Se incluye un script para la compilacion de OSQP y los codigos que implementan MPC, ademas se realiza la ejecucion de un ejemplo con un horizonte de prediccion N=4. Para ejecutar este script se deben otorgar permisos de ejecucion por medio de:

  `chmod +x setup_osqp.sh`
  
Luego, ejecutar el script por medio de:

  `./setup_osqp.sh`

Este script compilara los archivos del solver, eliminara algunas carpetas para evitar conflictos entre los makefile y finalmente ejecutara el archivo makefile del codigo MPC a modo de ejemplo. Si todo ha salido bien, la compilacion de archivos del solver y ejecucion de los codigos MPC se puede visualizar por medio de la terminal. Una salida de referencia generada por los codigos MPC es la siguiente:

**salida de referencia**

## Modo de uso - Ejecucion, muestras(samples) y salida de datos

Una vez ejecutado el script y se ha verificado la compilacion y ejecucion del ejemplo por terminal, se obtendra el ejecutable **runner** con horizonte de prediccion N=4 por defecto. Para ejecutar este programa con distintas configuraciones de muestras, salidas de datos e iteraciones, se debe utilizar el siguiente comando:

`./runner samplesMPC_N[HOR_PREDIC_MUESTRAS].txt [OUTPUT_DATOS] [OUTPUT_TIEMPO] [NUM_ITERACIONES] [NUMERO_EJECUCION]`

Donde:

samplesMPC_N[HOR_PREDIC_MUESTRAS].txt - Es el archivo de muestras obtenidas de MATLAB, HOR_PREDIC_MUESTRAS indica el horizonte de prediccion utilizado para crear las muestras

[OUTPUT_DATOS]: Valor binario para habilitar(1) o deshabilitar(0) la salida de datos como actuacion y variables de estado del sistema en un archivo csv.

[OUTPUT_TIEMPO]: Valor binario para habilitar(1) o deshabilitar(0) la salida de datos como tiempo de ejecucion del solver y lazo MPC en un archivo csv.

[NUM_ITERACIONES] - Parametro para indicar el maximo numero de iteraciones permitidas

[NUMERO_EJECUCION] - Parametro que agrega un numero o tag a los resultados de salidas o tiempo de los codigos

En caso de ser necesario cambiar la precision entre single y double o cambiar el horizonte de prediccion, antes de compilar se debe modificar en el header file *source/specs.h* el horizonte de predicción que se utilizará `#define N_QP 4` y se debe descomentar la linea `#define ELEM_FLOAT` o `#define ELEM_DOUBLE` para fijar el tamaño de punto flotante que se va a ocupar. Para compilar nuevamente se debe correr el makefile por medio de:

 `make`
 
Y para ejecutar se utiliza el comando definido anteriormente.
 
No olvidar que cada ejecucion de codigo debe ir junto a una archivo de muestras correspondiente al horizonte de prediccion definido, samplesMPC_N4.txt en el caso de ejemplo.

### Archivos de salida

Los archivos de salida generados tienen por nombre csv_qpoasesMPC[N]x[N]OUT_iter[NUM_ITERACIONES]_[NUMERO_EJECUCION]
 csv_qpoasesMPC[N]x[N]TIME_iter[NUM_ITERACIONES]_[NUMERO_EJECUCION], nombre de archivo para variables de estado y tiempos de ejecucion respectivamente.
 
Donde: 

[N]: Al igual que HOR_PREDIC_MUESTRAS, indica el horizonte de prediccion utilizado.

[NUM_ITERACIONES] - Indicador de maximo numero de iteraciones permitidas utilizado para generar los resultados.

[NUMERO_EJECUCION] - Indicador o tag de resultados de salidas o tiempo de la ejecucion de los codigos.

El formato del archivo de salidas(OUT) es el siguiente:

actuacion_1, vel.ang_1, pos.ang_1
actuacion_2, vel.ang_2, pos.ang_2
actuacion_3, vel.ang_3, pos.ang_3

El formato del archivo de tiempo(TIME) es el siguiente:

tiempo_lazo_1, tiempo_solver_1
tiempo_lazo_2, tiempo_solver_2
tiempo_lazo_3, tiempo_solver_3

### Creacion de muestras(samples)

Se tienen que tener las muestras generadas por el código *servoMPCReferenceTracking.m* con el horizonte de predicción que se va a utilizar. Para generar estas muestras seguir las instrucciones en GenerateSamplesMATLAB/Readme.md

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