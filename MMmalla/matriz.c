#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Bloques máximo de particionamiento de la matriz
#define rmax 4

//Numero máximos de elementos contenidos en el bloque
#define maxbloqtam 100


// Obtiene el resultado de si es un cuadrado perfecto y almacena el valor en r
short cuadradoPerfecto(const int *numprocs, int *r){

	float res= sqrt(*numprocs);
	int resultado=(int) res;

	if(res-resultado)
		return -1;
	(*r)=res;
	return 0;
}
// Multplica dos vectores como si fuesen matrices
void mult(double a[], double b[], double c[], int m) {
	int i,j,k;
	for (i=0; i<m; i++)
		for (j=0; j<m; j++)
			for (k=0; k<m; k++)
			{
				
				c[i*m+j]=c[i*m+j]+a[i*m+k]*b[k*m+j];
			}
	return; 
}
// Rellena el vector de A del proceso N
void rellenarVector(double *B,double *A, double *C,const int *bloqtam, const int *fila, const int *columna){
	int longVector=(*bloqtam)*(*bloqtam);
	int mfila= *fila;
	int mcolumna=*columna;
	int index=0;

	for(index=0; index<longVector; ++index){
		A[index]=index*(float)((mfila*mcolumna+1)*(mfila*mcolumna+1))/longVector;
		C[index]=0;
		B[index]=0;
	
	}
}
// Calcula la matriz identidad del proceso N
void matrizIdentidad(double B[],const int *bloqtam)
{
	int esUno=0;
	int index=0;
	int longVector=(*bloqtam)*(*bloqtam);
	int numElemento= (*bloqtam)+1;
	for(index=0; index<longVector; ++index)
	{
	if(esUno==0){
		B[index]=1;
		esUno=numElemento;
	}
	else
		B[index]=0;

	--esUno;
	}
}
	

int main(int argc, char **argv)
{
	//////////////// DECLARAMOS VARIABLES /////////////////////////////
	int numprocs, myrank, bloqtam, fila,columna, r,sizeBuffer,totalSize;
	//int procEnvio,procsAux,rangoProc;
	int numErrorLocal=0,numErrorTotal=0;
	double *A,*B,*ATMP,*C,*buffer;
	short index=0;
	int falloDeInicio=0;
	MPI_Init(&argc,&argv);

	///////////////// INICIAMOS MPI////////////
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
 	MPI_Status status;

 	//////// TRATAMIENTO DE ARGUMENTOS ////////
 	if(myrank==0)
 	{
 		if(argc < 2)
 		{
 			printf("Debes indicar el tamanyo de los bloques de la matriz \n");
 			falloDeInicio=-1;
 		}
 		else
			bloqtam=atoi(argv[1]);
 		if(bloqtam<=0 || bloqtam>maxbloqtam)
		{
			printf("tamanyo de los bloques de la matriz no mayor de 100\n");
			falloDeInicio=-1;
		}
 		if(cuadradoPerfecto(&numprocs,&r)==-1 || r>rmax )
 		{
			printf("La cantidad de procesos debe ser igual a la de un cuadrado perfecto (4,9,16,25...)\n");
			falloDeInicio=-1;
 		}
 	}
 	MPI_Bcast(&falloDeInicio,1,MPI_INT,0,MPI_COMM_WORLD);
 	if(falloDeInicio!=-1)
 	{	

 		////INICIALIZAMOS bloqtam y r en todos los procesos ////
	 	MPI_Bcast(&bloqtam,1,MPI_INT,0,MPI_COMM_WORLD);
	 	MPI_Bcast(&r,1,MPI_INT,0,MPI_COMM_WORLD);

	 	/// DECLARAMOS VARIABLES A UTILIZAR POR TODOS LOS PROCESOS ////
	 	int i=0;
	 	int posEnFila=0;
		//Obtenemos las posiciones de los procesos que se encuentran arriba y abajo
		int arriba=(myrank-r+numprocs)%numprocs;
		int abajo=(myrank+r+numprocs)%numprocs;
		int mifila[r]; // Contienen r elementos

		/// LOCALIZAMOS FILA, COLUMNA Y CALCULAMOS EL TAMAÑO DEL VECTOR /////
	 	fila= myrank/r;
	 	columna= myrank%r;
	 	totalSize= bloqtam*bloqtam;
		mifila[columna]=myrank; // almacenamos en la posicion de la fila el proceso
	 	 

	 	//Obtenemos las posiciones en la malla de los procesos
		for(i=0; i<numprocs; i++)
			if(i!=myrank)
				MPI_Send(&columna,1,MPI_INT,i,fila,MPI_COMM_WORLD);

		for(i=0; i<r-1;++i){
			MPI_Recv(&posEnFila,1,MPI_INT,MPI_ANY_SOURCE,fila,MPI_COMM_WORLD,&status);
			mifila[posEnFila]= status.MPI_SOURCE;
		}

		////// RESERVA DE MEMORIA PARA LOS DISTINTOS VECTORES QUE ENVIAREMOS ////
	 	//Cada proceso reserva el espacio necesario de cada vector
	 	B=(double *)malloc(totalSize* sizeof(double));
	 	A=(double *)malloc(totalSize* sizeof(double));
	 	ATMP=(double *)malloc(totalSize* sizeof(double));
	 	C= (double *)malloc(totalSize* sizeof(double));

	 	////INIZIALIZACION DE LOS ELEMENTOS DE LA MATRIZ /////
	 	// Inicializa los valores de B y A;
	 	rellenarVector(B,A,C,&bloqtam, &fila, &columna);
	 	if(myrank%(r+1)==0)
	 		matrizIdentidad(B, &bloqtam);



		////RESERVAMOS ESPACIO PARA EL BUFFER //////
		//Preparamos el buffer de envio de la matriz
		MPI_Pack_size(totalSize,MPI_DOUBLE,MPI_COMM_WORLD, &sizeBuffer);
		sizeBuffer =numprocs *(sizeBuffer +MPI_BSEND_OVERHEAD );
		buffer= (double *)malloc(sizeBuffer);
		MPI_Buffer_attach(buffer,sizeBuffer);
	 	

	 	///// BUCLE PRINCIPAL DEL ALGORITMOS //////
	 	for(index=0; index<r; ++index)
	 	{
	 		if(columna==((fila+index)%r))
	 		{
	 			for(i=0; i<r;++i)
	 				if(mifila[i]!=myrank)
	 					MPI_Bsend(A,totalSize,MPI_DOUBLE,mifila[i],mifila[i],MPI_COMM_WORLD);
	 			
	 			mult(A,B,C,bloqtam);
	 		}
	 		else
	 		{
				MPI_Recv(ATMP,totalSize,MPI_DOUBLE,mifila[(fila+index)%r],myrank,MPI_COMM_WORLD,&status);
				mult(ATMP,B,C,bloqtam);
	 		}
	 		MPI_Bsend(B,totalSize,MPI_DOUBLE,arriba,9,MPI_COMM_WORLD);
	 		MPI_Recv(B,totalSize,MPI_DOUBLE,abajo,9,MPI_COMM_WORLD,&status);	
	 	}
	 	
	 	///// CALCULO DE LOS ERRORES COMETIDO //////
	 	for(i=0; i<totalSize; ++i)
	 		if(A[i]!=C[i])
	 			++numErrorLocal;

	 	//// ACULUMAMOS EL CALCULO DEL ERROR EN EL PROCESO 0 PARA MOSTRAR RESULTADO FINAL /////
		MPI_Reduce(&numErrorLocal,&numErrorTotal,1,MPI_DOUBLE,MPI_SUM,0, MPI_COMM_WORLD);

		if(myrank==0)
			printf("NUMERO DE ERRORES : %d \n",numErrorTotal );


		/// LIBERAMOS MEMORIA RESERVADA DURANTE EL PROGRAMA /////
		MPI_Buffer_detach(buffer, &sizeBuffer);	
	 	free(A);
	 	free(B);
	 	free(C);
	 	free(ATMP);
	 	free(buffer);	
 	}


 	MPI_Finalize();
 	return 0;	

}


	