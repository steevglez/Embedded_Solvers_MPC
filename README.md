# Embedded_Solvers_MPC

#     Implementación de Algoritmo de Control Predictivo para Motor de Corriente Continua Utilizando Sistemas Embebidos


En este repositorio se encuentran los códigos asociados al trabajo de memoria "Implementación de Algoritmo de Control Predictivo para Motor de Corriente Continua Utilizando Sistemas Embebidos".

Se tiene código en MATLAB que se utiliza como referencia del correcto funcionamiento de MPC, se encuentra en la carpeta [GenerateSamplesMATLAB](GenerateSamplesMATLAB).
Es una version que incluye una modificacion que almacena las matrices A y B y del sistema, la version original puede encontrarse en https://github.com/morrisort/embeddedMPC, por Andrew Morrison. Ademas, se incluye la carpeta [Samples](Samples) que contiene las muestras de datos utilizadas en este trabajo, considerando horizontes de prediccion N=[2,3,4,5,6].

Los codigos confeccionados en este trabajo se encuentran en las carpetas [qpOASES-FIXED_ITER-N4](qpOASES-FIXED_ITER-N4), [OSQP-FIXED_ITER-N4](OSQP-FIXED_ITER-N4), [FiOrdOS-FIXED_ITER-N4](FiOrdOS-FIXED_ITER-N4) que contienen los codigos que implementan un lazo MPC en modelo de motor DC e integran el solver indicado. Cada carpeta esta asociada a un solver con la version mas optimizada lograda con un horizonte de prediccion N=4 a modo de ejemplo. 

---

## Replicar resultados

Para replicar los resultados se debe 

### Código Matlab



---
