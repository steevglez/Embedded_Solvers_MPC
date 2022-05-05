# Embedded_Solvers_MPC

#     Implementación de Algoritmo de Control Predictivo para Motor de Corriente Continua Utilizando Sistemas Embebidos


En este repositorio se encuentran los códigos asociados al trabajo de memoria "Implementación de Algoritmo de Control Predictivo para Motor de Corriente Continua Utilizando Sistemas Embebidos".

La carpeta [Matlab](Matlab) contiene archivos de Matlab de referencia en la carpeta [MATLAB_MPC_ref](Matlab/MATLAB_MPC_ref) utilizado como referencia del correcto funcionamiento de MPC. Tambien, se incluyen los archivos utilizados para crear las muestras de datos en [GenerateSamplesMATLAB](Matlab/GenerateSamplesMATLAB), esta una version donde se incluye una modificacion que almacena las matrices A y B del sistema fisico, la version original puede encontrarse en https://github.com/morrisort/embeddedMPC, por Andrew Morrison. Ademas, se incluye la carpeta [Samples](Matlab/samples) que contiene las muestras de datos utilizadas en este trabajo, considerando horizontes de prediccion N=[2,3,4,5,6].

Los codigos confeccionados en este trabajo se encuentran en las carpetas [qpOASES-FIXED_ITER-N4](qpOASES-FIXED-ITER-N4), [OSQP-FIXED_ITER-N4](OSQP-FIXED-ITER-N4), [FiOrdOS-FIXED_ITER-N4](FiOrdOS-FIXED-ITER-N4) que contienen los codigos que implementan un lazo MPC en modelo de motor DC e integran el solver indicado. Cada carpeta esta asociada a un solver con la version mas optimizada lograda con un horizonte de prediccion N=4 a modo de ejemplo. 

---

## Replicar resultados

Para replicar los resultados se debe 

### Código Matlab



---
