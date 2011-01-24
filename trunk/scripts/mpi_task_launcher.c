#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main( argc, argv )
int  argc;
char **argv;
{
    int rank, size;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    rank++;

    printf( "Launching task %d of %d\n", rank, size );

    char cmd[2048];
    sprintf(cmd,"TASK_ID=%i %s > out.p%i 2>&1",rank,argv[1],rank);
    system(cmd);

    MPI_Finalize();
    return 0;
}
