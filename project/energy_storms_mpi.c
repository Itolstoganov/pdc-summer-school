/*
 * Simplified simulation of high-energy particle storms
 *
 * Parallel computing (Degree in Computer Engineering)
 * 2017/2018
 *
 * Version: 2.0
 *
 * Code prepared to be used with the Tablon on-line judge.
 * The current Parallel Computing course includes contests using:
 * OpenMP, MPI, and CUDA.
 *
 * (c) 2018 Arturo Gonzalez-Escribano, Eduardo Rodriguez-Gutiez
 * Grupo Trasgo, Universidad de Valladolid (Spain)
 *
 * This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.
 * https://creativecommons.org/licenses/by-sa/4.0/
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

/* Headers for the MPI assignment versions */
#include<mpi.h>

/* Use fopen function in local tests. The Tablon online judge software 
   substitutes it by a different function to run in its sandbox */
#ifdef CP_TABLON
#include "cputilstablon.h"
#else
#define    cp_open_file(name) fopen(name,"r")
#endif

/* Function to get wall time */
double cp_Wtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}


#define THRESHOLD    0.001f

/* Structure used to store data for one storm of particles */
typedef struct {
    int size;    // Number of particles
    int *posval; // Positions and values
} Storm;

/* THIS FUNCTION CAN BE MODIFIED */
/* Function to update a single position of the layer */
void update( float *layer, int layer_size, int k, int k_offset, int pos, float energy ) {
    /* 1. Compute the absolute value of the distance between the
        impact position and the k-th position of the layer */
    int distance = pos - (k + k_offset);
    if ( distance < 0 ) distance = - distance;

    /* 2. Impact cell has a distance value of 1 */
    distance = distance + 1;

    /* 3. Square root of the distance */
    /* NOTE: Real world atenuation typically depends on the square of the distance.
       We use here a tailored equation that affects a much wider range of cells */
    float atenuacion = sqrtf( (float)distance );

    /* 4. Compute attenuated energy */
    float energy_k = energy / layer_size / atenuacion;

    /* 5. Do not add if its absolute value is lower than the threshold */
    if ( energy_k >= THRESHOLD / layer_size || energy_k <= -THRESHOLD / layer_size )
        layer[k] = layer[k] + energy_k;
}


/* ANCILLARY FUNCTIONS: These are not called from the code section which is measured, leave untouched */
/* DEBUG function: Prints the layer status */
void debug_print(int layer_size, float *layer, int *positions, float *maximum, int num_storms ) {
    int i,k;
    /* Only print for array size up to 35 (change it for bigger sizes if needed) */
    if ( layer_size <= 35 ) {
        /* Traverse layer */
        for( k=0; k<layer_size; k++ ) {
            /* Print the energy value of the current cell */
            printf("%10.4f |", layer[k] );

            /* Compute the number of characters. 
               This number is normalized, the maximum level is depicted with 60 characters */
            int ticks = (int)( 60 * layer[k] / maximum[num_storms-1] );

            /* Print all characters except the last one */
            for (i=0; i<ticks-1; i++ ) printf("o");

            /* If the cell is a local maximum print a special trailing character */
            if ( k>0 && k<layer_size-1 && layer[k] > layer[k-1] && layer[k] > layer[k+1] )
                printf("x");
            else
                printf("o");

            /* If the cell is the maximum of any storm, print the storm mark */
            for (i=0; i<num_storms; i++) 
                if ( positions[i] == k ) printf(" M%d", i );

            /* Line feed */
            printf("\n");
        }
    }
}

/*
 * Function: Read data of particle storms from a file
 */
Storm read_storm_file( char *fname ) {
    FILE *fstorm = cp_open_file( fname );
    if ( fstorm == NULL ) {
        fprintf(stderr,"Error: Opening storm file %s\n", fname );
        exit( EXIT_FAILURE );
    }

    Storm storm;    
    int ok = fscanf(fstorm, "%d", &(storm.size) );
    if ( ok != 1 ) {
        fprintf(stderr,"Error: Reading size of storm file %s\n", fname );
        exit( EXIT_FAILURE );
    }

    storm.posval = (int *)malloc( sizeof(int) * storm.size * 2 );
    if ( storm.posval == NULL ) {
        fprintf(stderr,"Error: Allocating memory for storm file %s, with size %d\n", fname, storm.size );
        exit( EXIT_FAILURE );
    }
    
    int elem;
    for ( elem=0; elem<storm.size; elem++ ) {
        ok = fscanf(fstorm, "%d %d\n", 
                    &(storm.posval[elem*2]),
                    &(storm.posval[elem*2+1]) );
        if ( ok != 2 ) {
            fprintf(stderr,"Error: Reading element %d in storm file %s\n", elem, fname );
            exit( EXIT_FAILURE );
        }
    }
    fclose( fstorm );

    return storm;
}

/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
    int i,j,k;

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* 1.1. Read arguments */
    if (argc<3) {
        fprintf(stderr,"Usage: %s <size> <storm_1_file> [ <storm_i_file> ] ... \n", argv[0] );
        exit( EXIT_FAILURE );
    }

    int layer_size = atoi( argv[1] );
    int num_storms = argc-2;
    Storm storms[ num_storms ];

    /* 1.2. Read storms information */
    for( i=2; i<argc; i++ ) {
        storms[i-2] = read_storm_file( argv[i] );
    }

    /* 1.3. Intialize maximum levels to zero */
    float maximum[ num_storms ];
    int positions[ num_storms ];
    for (i=0; i<num_storms; i++) {
        maximum[i] = 0.0f;
        positions[i] = 0;
    }

    /* 2. Begin time measurement */
    MPI_Barrier(MPI_COMM_WORLD);
    double ttotal = cp_Wtime();

    /* START: Do NOT optimize/parallelize the code of the main program above this point */

    /* 3. Allocate memory for the layer and initialize to zero */
    int num_chunks = size - 1;
    int pos_offset = 0;
    int chunk_size = layer_size / num_chunks;
    int max_chunk_size = 2 * chunk_size;

    float *layer = NULL;
    float *layer_copy = NULL;

    float *chunk = (float *)malloc( sizeof(float) * max_chunk_size);
    for( k=0; k<max_chunk_size; k++ ) chunk[k] = 0.0f;

    /* 4. Storms simulation */
    if (rank == 0) {
        layer = (float *)malloc( sizeof(float) * layer_size );
        layer_copy = (float *)malloc( sizeof(float) * layer_size );
        if ( layer == NULL || layer_copy == NULL ) {
            fprintf(stderr,"Error: Allocating the layer memory\n");
            exit( EXIT_FAILURE );
        }
        for( k=0; k<layer_size; k++ ) layer[k] = 0.0f;
        for( k=0; k<layer_size; k++ ) layer_copy[k] = 0.0f;

        pos_offset = 0;
        for (int rec_rank = 1; rec_rank <= num_chunks; ++rec_rank) {
            int current_chunk_size = rec_rank < num_chunks ? chunk_size : layer_size - pos_offset;
            MPI_Send(&current_chunk_size, 1, MPI_INT, rec_rank, 2 * rec_rank, MPI_COMM_WORLD);
            MPI_Send(&pos_offset, 1, MPI_INT, rec_rank, 2 * rec_rank + 1, MPI_COMM_WORLD);
            pos_offset += current_chunk_size;
        }
    } else {
        MPI_Status status;
        MPI_Recv(&chunk_size, 1, MPI_INT, 0, 2 * rank, MPI_COMM_WORLD, &status);
        MPI_Recv(&pos_offset, 1, MPI_INT, 0, 2 * rank + 1, MPI_COMM_WORLD, &status);
    }

    for( i=0; i<num_storms; i++) {
        if (rank == 0) {
            for( k=0; k<layer_size; k++ ) {
                layer_copy[k] = 0;
            }

            pos_offset = 0;
            for (int rec_rank = 1; rec_rank <= num_chunks; ++rec_rank) {
                int recieved;
                int current_chunk_size = rec_rank < num_chunks ? chunk_size : layer_size - pos_offset;
                MPI_Status status;
                MPI_Recv(layer_copy + pos_offset, current_chunk_size, MPI_FLOAT, rec_rank, (2 + i) * num_chunks + rec_rank, MPI_COMM_WORLD, &status);
                // MPI_Get_count(&status, MPI_INT, &recieved);
                // fprintf(stdout, "Storm %d, received %d\n", i, recieved);
                for (int k = pos_offset; k < pos_offset + current_chunk_size; ++k) {
                    // fprintf(stdout, "%f, ", layer_copy[k]);
                    layer[k] += layer_copy[k];
                }
                // fprintf(stdout, "\n");
                pos_offset += current_chunk_size;
            }

            // fprintf(stdout, "Storm %d, relaxation\n", i);
            /* 4.2. Energy relaxation between storms */
            /* 4.2.1. Copy values to the ancillary array */
            for( k=0; k<layer_size; k++ ) {
                layer_copy[k] = layer[k];
            }

            /* 4.2.2. Update layer using the ancillary values.
                      Skip updating the first and last positions */
            for( k=1; k<layer_size-1; k++ )
                layer[k] = ( layer_copy[k-1] + layer_copy[k] + layer_copy[k+1] ) / 3;

            /* 4.3. Locate the maximum value in the layer, and its position */
            for( k=1; k<layer_size-1; k++ ) {
                /* Check it only if it is a local maximum */
                if ( layer[k] > layer[k-1] && layer[k] > layer[k+1] ) {
                    if ( layer[k] > maximum[i] ) {
                        maximum[i] = layer[k];
                        positions[i] = k;
                    }
                }
            }
        } else {
            MPI_Status status;
            /* 4.1. Add impacts energies to layer cells */
            /* For each particle */
            for( j=0; j<storms[i].size; j++ ) {
                /* Get impact energy (expressed in thousandths) */
                float energy = (float)storms[i].posval[j*2+1] * 1000;
                /* Get impact position */
                int position = storms[i].posval[j*2];

                /* For each cell in the layer */
                for( k=0; k<chunk_size; k++ ) {
                    /* Update the energy value for the cell */
                    update( chunk, layer_size, k, pos_offset, position, energy );
                }
            }
            // fprintf(stdout, "Storm %d, sending from %d, position %d, chunk size %d\n", i, rank, pos_offset, chunk_size);
            MPI_Send(chunk, chunk_size, MPI_FLOAT, 0, (2 + i) * num_chunks + rank, MPI_COMM_WORLD);
            for (k = 0; k < chunk_size; ++k) {
                // fprintf(stdout, "%f, ", chunk[k]);
                chunk[k] = 0;
            }
            // fprintf(stdout, "\n");
        }
    }   

    /* END: Do NOT optimize/parallelize the code below this point */

    /* 5. End time measurement */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {

        ttotal = cp_Wtime() - ttotal;

        /* 6. DEBUG: Plot the result (only for layers up to 35 points) */
        #ifdef DEBUG
        debug_print( layer_size, layer, positions, maximum, num_storms );
        #endif

        /* 7. Results output, used by the Tablon online judge software */
        printf("\n");
        /* 7.1. Total computation time */
        printf("Time: %lf\n", ttotal );
        /* 7.2. Print the maximum levels */
        printf("Result:");
        for (i=0; i<num_storms; i++)
            printf(" %d %f", positions[i], maximum[i] );
        printf("\n");


        /* 8. Free resources */    
        for( i=0; i<argc-2; i++ )
            free( storms[i].posval );
    }

    /* 9. Program ended successfully */
    MPI_Finalize( );
    return 0;
}

