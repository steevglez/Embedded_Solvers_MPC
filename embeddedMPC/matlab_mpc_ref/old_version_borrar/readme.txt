Los archivos servo_mpc_reference_tracking_dense.m y servo_mpc_reference_tracking_sparse.m contienen ejemplos de MPC para un servomotor de referencia utilizando las formulaciones densas y sparse para el problema QP, respectivamente.

El archivo pdip.m contiene la implementacion del algoritmo de punto interior para resolver el problema QP. Dentro de este problema se requieren resolver sistemas de ecuaciones lineales, para lo cual se pueden utilizar distintos algoritmos.

Los archivos lschol.m, min_res.m, y cg_wip.m resuleven sistemas de ecuaciones lineales utilizando Cholesky, Minres, y conjungate gradient, respectivamente.

El objetivo del trabajo es implementar en FPGA alguna tecnica para resolver sistemas de ecuaciones lineales, que pueden ser las incluidas aca u otras.

Resolver sistemas de ecuaciones lineales es solo una de las tareas a resolver para implementar MPC. Para obtener sistemas de ecuaciones de ejemplo y soluciones obtenidas con el codigo en matlab para generar golden references para testbenches, se deben guardar las matrices Ak y bk que ingresan a la parte del solver en el codigo pdip.m . 


Notar que para un mismo problema, las formulaciones densas y sparse generaran matrices de distinto orden y caracteristicas para el linear solver. El objetivo es explorar los tradeoff asociadas a cada formulacion del problema al momento de implementar el linear solver en hardware.

