#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>


typedef struct _MinMax{
   double Min;
   double Max;
} MinMax;

int main(int argc, char **argv)
{
	int **Crear_matriz(int , int );
	void Liberar_matriz(int **, int );
	MinMax BuscarMinMax(int **, int , int);
	double ctimer(void);

	int pixelX,pixelY;               // Coordenadas de la imagen
	int pixelXmax = 1024; 
	int pixelYmax = 1024;
	double Creal,Cimg;               // Coordenadas de los puntos complejos
	const double RealMinArray[] = {-2, -1.023438, -1.017579, -1.017523, -1.0190739863281251, -1.0184326403503419, -1.0175024721407624, -1.0176950869990224, 0.2720};
	const double RealMaxArray[] = { 1, -0.992188, -1.016968, -1.017493, -1.0178532832031251, -1.0184266798858643, -1.0175010058651026, -1.0173896129032258, 0.3720};
	const double ImMinArray[] = {-1.5, -0.285156, -0.274444, -0.274065, -0.2672003476562500, -0.2667537630310058,  -0.2740544516129032, -0.2772175483870968, 0.4805};
	const double ImMaxArray[] = { 1.5, -0.25,     -0.273758, -0.274032, -0.2658270664062500, -0.2667336466064453,  -0.2740528387096774, -0.2768738924731183, 0.5805}; 

	double RealMin=-2.0; 
	double RealMax=1.0;  
	double ImMin=-1.5; 
	double ImMax=1.5; 
	double AnchoPixel; 
	double AltoPixel;
	double Tinicial, Tfinal, Ttotal;

	const int MaxValorTonos=255; 
	FILE *ImgFile, *ImgFile2;
	char ArchivoImagen[]="12345678901234567890imgA.pgm";
	char ArchivoImagen2[]="12345678901234567890imgB.pgm";
	char *comentario="# ";
	int bn, bn2;
	double Zx, Zy;             // Z=Zx+Zy*i  Zx parte real, Zy parte imaginaria
	double Zx2, Zy2;           // Zx2=Zx*Zx, Zy2=Zy*Zy  
	int Iter,i,j;
	int IterMax=1000;
	const double Salida=2;    // valor de escape
	double Salida2=Salida*Salida, SumaExponencial;
	int dominio=0;
	int **matriz, **matriz2;
	MinMax  min_max_img;
	int myrank, numprocs;
	int fila;
	MPI_Status estado, status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


	//printf("Soy el proceso %d de un total de %d.\n",myrank,numprocs);
	if (myrank == 0)
	{
		printf("\n  *************** POSIBLES PARÁMETROS ***************\n");
		printf("  * Ancho imagen: pixelXmax                         *\n");
		printf("  * Alto imagen: pixelYmax                          *\n");
		printf("  * Numero de iteraciones: IterMax                  *\n");
		printf("  * Nombre archivo imagen de salida: ArchivoImagen  *\n");
		printf("  * Intervalo a considerar: dominio                 *\n");  
		printf("  ***************************************************\n\n");
		switch (argc)
		{
		    case 6: sscanf(argv[5], "%i", &dominio);
		            if (8*(dominio + 1) > sizeof(RealMinArray)){
		               printf("El array de dominios no es tan grande\n");
		               return 0;
		            }
		            RealMin = RealMinArray[dominio];
		            RealMax = RealMaxArray[dominio];
		            ImMin = ImMinArray[dominio];
		            ImMax = ImMaxArray[dominio];
		    case 5: strcpy(ArchivoImagen, argv[4]);
		            strcat(ArchivoImagen, "A.pgm");
		            strcpy(ArchivoImagen2, argv[4]);
		            strcat(ArchivoImagen2, "B.pgm");
		    case 4: sscanf(argv[3], "%i", &IterMax);
		    case 3: sscanf(argv[2], "%i", &pixelYmax);
		    case 2: sscanf(argv[1], "%i", &pixelXmax);
		    case 1: break;
		    default: printf("Demasiados parametros\n");
		             return 0;
		}

		printf("  *************** DATOS DE LA EJECUCIÓN ****************************\n");
		printf("  * pixelXmax = %4d, pixelYmax = %4d                             *\n", pixelXmax, pixelYmax);
		printf("  * IterMax = %5d                                                *\n", IterMax);
		printf("  * ArchivoImagen = %15s %15s                *\n", ArchivoImagen, ArchivoImagen2);
		printf("  * Dominio = %2d                                                   *\n", dominio);
		printf("  * Intervalo para X: [%20.15f,%20.15f]  *\n  * Intervalo para Y: [%20.15f,%20.15f]  *\n", RealMin, RealMax, ImMin, ImMax);
		printf("  ******************************************************************\n\n");
		/* Se crea un nuevo archivo y se abre en binario  */
		ImgFile= fopen(ArchivoImagen,"wb"); 
		ImgFile2= fopen(ArchivoImagen2,"wb"); 
		/* Se escribe la cabecera - ASCII  */
		fprintf(ImgFile,"P5\n %s\n %d\n %d\n %d\n",comentario,pixelXmax,pixelYmax,MaxValorTonos);
		fprintf(ImgFile2,"P5\n %s\n %d\n %d\n %d\n",comentario,pixelXmax,pixelYmax,MaxValorTonos);
	}

	MPI_Bcast(&pixelYmax,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&pixelXmax,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&IterMax,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&RealMin,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&RealMax,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ImMin,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ImMax,1,MPI_INT,0,MPI_COMM_WORLD);

	AnchoPixel=(RealMax-RealMin)/(pixelXmax-1); 
	AltoPixel=(ImMax-ImMin)/(pixelYmax-1);

	if(myrank==0)
	{
		matriz = Crear_matriz(pixelYmax, pixelXmax);
		matriz2 = Crear_matriz(pixelYmax, pixelXmax);
		Tinicial = MPI_Wtime();
		int k; /// almacenamos al fila que ha procesado el proceso N
		//// COMIENZA A ENVIAR FILAS DE LA MATRIZ A LOS PROCESOS ////
		for(fila=0; fila<numprocs-1;fila++){
			// printf("Soy proceso %d y envio a %d\n",myrank,fila+1);
			MPI_Send(&fila, 1, MPI_INT, fila+1, 1, MPI_COMM_WORLD);   
			//envia al proceso fila+1
			MPI_Send(&matriz[fila][0], pixelXmax, MPI_INT, fila+1, 0, MPI_COMM_WORLD);
			MPI_Send(&matriz2[fila][0], pixelXmax, MPI_INT, fila+1, 0, MPI_COMM_WORLD);
      	}

      	/// RECIBO DE CADA PROCESO LOS RESULTADO DE SUS CALCULOS Y LOS ALMACENO EN LA MATRIZ
      	/// POSTERIORMENTE ENVIO LA SIGUIENTES FILAS A CALCULAR
      	for(fila=numprocs-1;fila<pixelYmax; ++fila)
		{
			MPI_Recv(&k, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &estado);
			MPI_Recv(&matriz[k][0], pixelXmax, MPI_INT, estado.MPI_SOURCE, estado.MPI_SOURCE, MPI_COMM_WORLD, &status);
			MPI_Recv(&matriz2[k][0], pixelXmax, MPI_INT, estado.MPI_SOURCE, estado.MPI_SOURCE, MPI_COMM_WORLD, &status);
			// printf("Soy proceso %d y acabo de recibir el trabajo y envio a %d \n",myrank, status.MPI_SOURCE );
			MPI_Send(&fila, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);   
			//envia al proceso fila+1
			MPI_Send(&matriz[fila][0], pixelXmax, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			MPI_Send(&matriz2[fila][0], pixelXmax, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}

		//// FINALMENTE RECIBO LOS RESULTADO DE LOS PROCESOS QUE QUEDAN //// 
		for (i=1; i<numprocs; i++) {
         	MPI_Recv(&k, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &estado);
			MPI_Recv(&matriz[k][0], pixelXmax, MPI_INT, estado.MPI_SOURCE, estado.MPI_SOURCE, MPI_COMM_WORLD, &status);
			MPI_Recv(&matriz2[k][0], pixelXmax, MPI_INT, estado.MPI_SOURCE, estado.MPI_SOURCE, MPI_COMM_WORLD, &status);

      	}
      	//////MANDA FINALIZACION DE EJECUCÓN //////
      	// printf("Soy proceso %d y envio para que terminen\n",myrank );
		for (i=1; i<numprocs; i++) {
			MPI_Send(&i, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
		}
		
	}
	else
	{
		//// EJECUTAMOS BUCLE INFINITO RECIBIENDO LAS FILAS QUE NOS MANDA 0 /////
		int *filaA,*filaB; 
		filaA = (int *)malloc(sizeof(int) * pixelXmax);
		filaB = (int *)malloc(sizeof(int) * pixelXmax);
		while (1) {
			MPI_Recv(&fila, 1, MPI_INT, 0,MPI_ANY_TAG, MPI_COMM_WORLD, &estado); 
			if (estado.MPI_TAG == 99) 
			{
				free(filaA);
				free(filaB);
				MPI_Finalize();
				return 0;
			}
			MPI_Recv(filaA, pixelXmax, MPI_INT, 0, 0, MPI_COMM_WORLD, &estado);
			MPI_Recv(filaB, pixelXmax, MPI_INT, 0, 0, MPI_COMM_WORLD, &estado);

			// IMPLEMENTACION DEL ALGORITMO PARA CALCULAR
			Cimg=ImMin + fila*AltoPixel;
			for(pixelX=0;pixelX<pixelXmax;pixelX++)
			{         
			   Creal=RealMin + pixelX*AnchoPixel;
			   Zx=0.0;         // Valor inicial
			   Zy=0.0;
			   Zx2=Zx*Zx;
			   Zy2=Zy*Zy;
			   SumaExponencial=0;
			   for (Iter=0;Iter<IterMax && ((Zx2+Zy2)<Salida2);Iter++)
			   {
			       SumaExponencial += exp( -sqrt(Zx2+Zy2) )/IterMax;   // se mantiene siempre entre (0,1)            
			         Zy=2*Zx*Zy + Cimg;     //parte imaginario
			       Zx=Zx2-Zy2 + Creal;    //parte real
			       Zx2=Zx*Zx;
			       Zy2=Zy*Zy;
			   }
			   if (Iter==IterMax) { /*  interior del conjunto Mandelbrot = negro */    
			      bn=0;          
			      bn2=0;            
			   }
			   else { /* exterior del conjunto Mandelbrot = blanco modificado con exp */
			      bn = MaxValorTonos - SumaExponencial * 255; 
			      bn2 = (Iter +1 - (int)(log(log(sqrt(Zx2+Zy2)) / log(2)) / log(2)))*255;
			   }
			   filaA[pixelX]=bn;
			   filaB[pixelX]=bn2;
			}
			

			// UNA VEZ FINALIZADO EL CALCULO ENVIO AL PROCESO 0 LOS RESULTADO
			// printf("Soy %d y envio mi info a 0\n", myrank);
			MPI_Send(&fila, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
			MPI_Send(filaA, pixelXmax, MPI_INT, 0, myrank, MPI_COMM_WORLD);
			MPI_Send(filaB, pixelXmax, MPI_INT, 0, myrank, MPI_COMM_WORLD);
		}
		
	}

	//////-----ALGORITMO DE PROCESAMIENTO DEL CONJUNTO -----_/////

	// for(pixelY=0;pixelY<pixelYmax;pixelY++)
	// {
	// 	Cimg=ImMin + pixelY*AltoPixel;
	// 	for(pixelX=0;pixelX<pixelXmax;pixelX++)
	// 	{         
	// 	   Creal=RealMin + pixelX*AnchoPixel;
	// 	   Zx=0.0;         // Valor inicial
	// 	   Zy=0.0;
	// 	   Zx2=Zx*Zx;
	// 	   Zy2=Zy*Zy;
	// 	   SumaExponencial=0;
	// 	   for (Iter=0;Iter<IterMax && ((Zx2+Zy2)<Salida2);Iter++)
	// 	   {
	// 	       SumaExponencial += exp( -sqrt(Zx2+Zy2) )/IterMax;   // se mantiene siempre entre (0,1)            
	// 	         Zy=2*Zx*Zy + Cimg;     //parte imaginario
	// 	       Zx=Zx2-Zy2 + Creal;    //parte real
	// 	       Zx2=Zx*Zx;
	// 	       Zy2=Zy*Zy;
	// 	   }
	// 	   if (Iter==IterMax) { /*  interior del conjunto Mandelbrot = negro */    
	// 	      bn=0;          
	// 	      bn2=0;            
	// 	   }
	// 	   else { /* exterior del conjunto Mandelbrot = blanco modificado con exp */
	// 	      bn = MaxValorTonos - SumaExponencial * 255; 
	// 	      bn2 = (Iter +1 - (int)(log(log(sqrt(Zx2+Zy2)) / log(2)) / log(2)))*255;
	// 	   }
	// 	   matriz[pixelY][pixelX] = bn;
	// 	   matriz2[pixelY][pixelX] = bn2;
	// 	}
	// }

	//////-----FIN ALGORITMO DE PROCESAMIENTO DEL CONJUNTO -----_/////
	if(myrank==0)
	{
		Tfinal = MPI_Wtime();
		Ttotal = Tfinal-Tinicial;
		printf("\nTiempo: %f segundos\n", Ttotal);
		min_max_img = BuscarMinMax(matriz2,pixelYmax,pixelXmax);
		for (i=0;i<pixelYmax;i++) {
		    for (j=0;j<pixelXmax;j++) {
		        matriz2[i][j] = matriz2[i][j] - min_max_img.Min;
		        matriz2[i][j] = matriz2[i][j] * (255.0/(min_max_img.Max-min_max_img.Min));
		    }
		}
		min_max_img = BuscarMinMax(matriz2,pixelYmax,pixelXmax);

		for (pixelY=0;pixelY<pixelYmax;pixelY++) {
		    for(pixelX=0;pixelX<pixelXmax;pixelX++){
		       fwrite(&matriz[pixelY][pixelX],1,1,ImgFile);
		       fwrite(&matriz2[pixelY][pixelX],1,1,ImgFile2);
		    }
		} 
		fclose(ImgFile);
		fclose(ImgFile2);
		Liberar_matriz(matriz,pixelYmax);
		Liberar_matriz(matriz2,pixelYmax);		
	}

	MPI_Finalize();
	return 0;
	
 }

int **Crear_matriz(int fila, int col)
{
    int **ret_val;
    int i;
 
    ret_val = (int **)malloc(sizeof(int *) * fila);
    if (ret_val == NULL) {
        perror("Problemas al dimensionar");
        exit(EXIT_FAILURE);
    }
 
    for (i = 0; i < fila; ++i) {
        ret_val[i] = (int *)malloc(sizeof(int) * col);
        if (ret_val[i] == NULL) {
            perror("Problemas al dimensionar");
            exit(EXIT_FAILURE);
        }
    }
    return ret_val;
}

void Liberar_matriz(int **matriz, int fila)
{
    int i;
 
    for (i = 0; i < fila; ++i)
        free(matriz[i]);
    free(matriz);
}

#include <sys/time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

static int nclock;
double ctimer(void)
{
   struct timeval tp;
   struct timezone tzp;
   double diff;
   nclock=sysconf(_SC_CLK_TCK);
   gettimeofday(&tp, &tzp);
   diff=(double)tp.tv_sec+(double)tp.tv_usec/1.0e6;
   return diff;
}

MinMax BuscarMinMax(int **array, int fil, int col) {
   MinMax  min_max;
   int index,fila,columna,fila_next,columna_next;
   int n = fil*col; 
   if ( n%2 != 0 ){

     min_max.Min = array[0][0];
     min_max.Max = array[0][0];

     index = 1;
   }
   else{
     if ( array[0][0] < array[0][1] ){
       min_max.Min = array[0][0];
       min_max.Max = array[0][1];
      }
      else {
       min_max.Min = array[0][1];
       min_max.Max = array[0][0];
      }
       index = 2;
   }

   int big, small,i;
   for ( i = index; i < n-1; i = i+2 ){
      fila = i / col;
      columna = i % col;  
      fila_next = fila;
      columna_next = columna + 1;
      if (columna_next == col) {
         fila_next = fila +1;
         columna_next = 0;
      }
      if ( array[fila][columna] < array[fila_next][columna_next] ){
        small = array[fila][columna];
        big = array[fila_next][columna_next];
      }
      else{
        small = array[fila_next][columna_next];
        big = array[fila][columna];
      }
      if ( min_max.Min > small ){
        min_max.Min = small;
      }
      if ( min_max.Max < big ){ 
        min_max.Max = big;
      }
   }
   printf("Minimo = %f, Maximo = %f\n", min_max.Min, min_max.Max);
   return min_max;
}










