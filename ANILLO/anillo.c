#include <stdio.h>
#include <mpi.h>


main (int argc, char **argv)
{
	int myrank, numprocs, resultlen,numero,i,nextProc;
	MPI_Status estado;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
 	nextProc=(myrank+1)%numprocs;

   	if(myrank==0)
   	{
   		double start_time, end_time;
   		printf("Soy proceso root, indica el numero:");
      scanf("%d", &numero);
   		start_time = MPI_Wtime();
      MPI_Send(&numero,1,MPI_INT,1,0,MPI_COMM_WORLD);
      
      MPI_Recv(&numero,1,MPI_INT,numprocs-1,0,MPI_COMM_WORLD,&estado);
      numero+=1;
      end_time= MPI_Wtime();
      printf("Tiempo en ejecutar:%f La suma total es: %d\n",end_time-start_time,numero );
   	}
   	else if(myrank!=0)
   	{
	   		MPI_Recv(&numero,1,MPI_INT,myrank-1,0,MPI_COMM_WORLD,&estado);
	   		numero+=1;
   			MPI_Send(&numero,1,MPI_INT,nextProc,0,MPI_COMM_WORLD);	
	}
	MPI_Finalize();
}