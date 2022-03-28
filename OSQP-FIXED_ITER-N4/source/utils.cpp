#include "utils.h"

void displayMatrix(elem* C, int rows, int cols){
	for (int r=0; r<rows; r++){
		for (int c=0; c<cols; c++){
			printf("%f ", C[r*cols + c]);
		}
		printf("\n");
	}
}

void genRandArray(double min, double max, int size, double *array){
	// generates random array. Values from min to max

	for(int i=0; i<size; i++){
		array[i] = min + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(max-min)));
	}
}

void genSymPosMat(double max, int size, double *matrix){
	// generate a random, symmetric and positive defined matrix
	// values from 0 to max
	double randNumber = 0;
	for (int r=0; r<size; r++){
		for(int c=r; c<size; c++){
			randNumber = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(max)));
			if (r == c){
				matrix[r*size + c] = randNumber;
			}
			else {
				matrix[r*size + c] = randNumber;
				matrix[c*size + r] = randNumber;
			}
		}
	}
}

int compare(double* gold, elem* result, int size, double th){
	int errors = 0;
	double dif = 0;
	for (int i=0; i<size; i++){
		dif = fabs(gold[i] - result[i]);
		// a comparison with NaN will always be false
		if (!(dif < th)){
			errors++;
		}
	}
	return errors;
}

int compare(float* gold, elem* result, int size, double th){
	int errors = 0;
	float dif = 0;
	for (int i=0; i<size; i++){
		dif = fabs(gold[i] - result[i]);
		// a comparison with NaN will always be false
		if (!(dif < th)){
			errors++;
		}
	}
	return errors;
}

void cast2elem(double *inputArray, elem *outputArray, int size){
	for(int i=0; i<size; i++){
		outputArray[i] = static_cast<elem>(inputArray[i]);
	}
}

void cast2float(double *inputArray, float *outputArray, int size){
	for(int i=0; i<size; i++){
		outputArray[i] = static_cast<float>(inputArray[i]);
	}
}

