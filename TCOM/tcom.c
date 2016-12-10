#include <stdio.h>
#include <mpi.h>


#define tam 15
#define microSec 1000000

main (int argc, char **argv)
{

	int myrank, 	//rango del proceso
	numprocs,		//numero de procesos
	i,j,			//iteradores
	sendProc,		//enviar al proceso siguiente o recibe del proceso anterior
	//valor en base a 2^N-> obtenemos el numero de elementos a enviar
	pot=0,			//numero de elementos a enviar 
	numIteraciones;	//Cantidad de iteraciones que se quieren realizar
	
	unsigned int tamanyos[tam]={256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304};
	unsigned int tamVector=4194304/sizeof(double);
	double datosEnvio[tamVector];
	MPI_Status estado;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
 	float beta,tau;
 	sendProc=(myrank+1)%numprocs;
 	//dato a enviar para comprobar latencia
 	char datoVacio;
 	//medida de tiempo
 	double start_time, end_time, time_total;

 	if(argc>=2)
 	{
 		
		sscanf(argv[1], "%i", &numIteraciones);
		if(numIteraciones<=0)
		{
			if(myrank==0)
				printf("Finaliza ejecución, numero de iteraciones menor o igual que 0\n");
			MPI_Finalize();
			return 0;
		}
 	}else
 	{
 		if(myrank==0)
			printf("Finaliza ejecución, indica el numero de iteraciones \n");
		MPI_Finalize();
		return 0;
 	}
 		



 	/////////////////////////////////////////////////////////
 	//medimos el tiempo que tarda en comunicarse dos procesos sin datos de transferencia
 	// obtenemos el parametro B
 	/////////////////////////////////////////////////////////
 	for(i=0; i<numIteraciones;++i)
 	{
 		if(myrank==0){
 
 		 	start_time = MPI_Wtime()*microSec;
 		 	// envia un mensaje al proces siguiente
 			MPI_Send(&datoVacio,1,MPI_BYTE,myrank+1,0,MPI_COMM_WORLD);
 			// espera un mensaje del proceso al que le ha enviado previamente el mensaje
 			MPI_Recv(&datoVacio,1,MPI_BYTE,myrank+1,0,MPI_COMM_WORLD,&estado);
 			end_time= MPI_Wtime()*microSec;
	 		time_total+=(end_time-start_time)/2.0;
			
 		}
 		else
 		{
 			//espera el mensaje del proceso anterior
 		 	MPI_Recv(&datoVacio,1,MPI_BYTE,myrank-1,0,MPI_COMM_WORLD,&estado);
 		 	// envia el mensaje al proceso del que le ha venido el mensaje
 			MPI_Send(&datoVacio,1,MPI_BYTE,myrank-1,0,MPI_COMM_WORLD);	
 		}
	}

	if(myrank==0){

		beta= time_total/=numIteraciones;
		printf("Tiempo de latencia en red (Beta): %f (microSec) \n",beta );
	}

	//Calculamos el tiempo de envio del array de doubles	
	/////////////////////////////////////////
	//	T= (tiempoEnvio - latencia )/ tamanyo 
	/////////////////////////////////////////
	for(i=0; i<tam; ++i){
			pot=tamanyos[i]/sizeof(double);

			for(j=0; j<numIteraciones; ++j){
				if(myrank==0){
		 		 	start_time = MPI_Wtime()*microSec;
		 		 	// envia un mensaje al proces siguiente
		 			MPI_Send(&datosEnvio,pot,MPI_DOUBLE,myrank+1,0,MPI_COMM_WORLD);
		 			// espera un mensaje del proceso al que le ha enviado previamente el mensaje
		 			MPI_Recv(&datosEnvio,pot,MPI_DOUBLE,myrank+1,0,MPI_COMM_WORLD,&estado);
		 			end_time= MPI_Wtime()*microSec;
			 		time_total+=(end_time-start_time)/2.0;
	 			}
	 			else
	 			{
	 				//espera el mensaje del proceso anterior
		 		 	MPI_Recv(&datosEnvio,pot,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD,&estado);
		 		 	// envia el mensaje al proceso del que le ha venido el mensaje
		 			MPI_Send(&datosEnvio,pot,MPI_DOUBLE,myrank-1,0,MPI_COMM_WORLD);	
	 			}
			}
			if(myrank==0){
				time_total/=numIteraciones;
				tau=(time_total-beta)/tamanyos[i];
				printf("numElementos: %d \ttamañoMensaje: %d (bytes) \tTau: %f (microSec) \tTcom: %f \n",pot,tamanyos[i],tau,tau*tamanyos[i]+beta);
			}
	}
	MPI_Finalize();
	return 0;

}