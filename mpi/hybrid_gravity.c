#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>

int Read_file(char *File_name, int Count, int *stars);
void calculate_gravitational_force(int Count, int *stars, double gravitational_force_list[], int my_rank, int size);
double calculate_distance(int *stars, int main_star_index, int other_star_index);
double calculatePartial_2Norm(int Count, double gravitational_force_list[], int size);
void write_file(double gravitational_force_list[], int Count);

int main(int argc, char *argv[])
{

    int size, my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    char *File_name;
    int Star_count;
    int *stars;

    int checkFile; //if file name is incorrect program doesn't run
    double start;//Clock
    
    double *gravitational_force_list;
    

    double normOfStars;       //Real Norm holding inside that variable
    double partialNorm = 0.0; // Partial norm hold in that

    // get variables from input
    File_name = argv[1];
    Star_count = atoi(argv[2]);

    //Allocate stars array pointer and graviatational force array with malloc
    stars = (int *)malloc(Star_count * 4 * sizeof(int));                             //Star İnformation Array
    double *partialForceList = (double *)malloc(Star_count / size * sizeof(double)); //Partial Gravitational Force List

    if (my_rank == 0)
    {
        gravitational_force_list = (double *)malloc(Star_count * sizeof(double)); //Gravitational Force List will collect inside this

        //Get Stars From File
        checkFile = Read_file(File_name, Star_count, stars);

		start = clock();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(stars, Star_count * 4, MPI_INT, 0, MPI_COMM_WORLD); // stars are scattering to processes
    MPI_Bcast(&checkFile, 1, MPI_INT, 0, MPI_COMM_WORLD);         // if file is incorrect scattering information to other processes

    if (checkFile == -1)
    {
        return 0;
    }

    //Calculate Gravitational Force
    calculate_gravitational_force(Star_count, stars, partialForceList, my_rank, size);

    //Forces are collecting from other process's partial forces in Process0
    MPI_Gather(partialForceList, Star_count / size, MPI_DOUBLE, gravitational_force_list, Star_count / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Calculate Partial 2 Norm with Partial Forces
    partialNorm = calculatePartial_2Norm(Star_count, partialForceList, size);

    MPI_Reduce(&partialNorm, &normOfStars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //get partial norm from other processes

    if (my_rank == 0)
    {

        printf("Calculating Time is : %lf  Seconds\n", (clock() - start) / 1000000);

        printf("2-Norm :  %2.13E \n", sqrt(normOfStars));//Print 2Norm

        //Write to file
        write_file(gravitational_force_list, Star_count);
        free(gravitational_force_list);
    }

    free(stars);
    free(partialForceList);

    MPI_Finalize();

    return 0;
}

//Get Stars From File
int Read_file(char *File_name, int Count, int *stars)
{

    FILE *file;
    file = fopen(File_name, "r");

    if (file)
    {

        for (int i = 0; i < Count; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                fscanf(file, "%d", (stars + i * 4 + j));
            }
        }
        printf("Star informations were got from file ! \n\n");
        fclose(file);
        return 1;
    }
    else
    {
        printf("File can not be opened !! \n");
    }

    return -1;
}

//Calculate Fravitational force
void calculate_gravitational_force(int Count, int *stars, double gravitational_force_list[], int my_rank, int size)
{

    const double gravitational_constat = 0.0006674; //Gravitational Constant (G)                                

    //Filling the pointer inside
    for (int i = 0; i < Count / size; i++)
    {
        gravitational_force_list[i] = 0;
    }

    printf("Calculating From %d.Process ...  ( the process may take a long time ) \n", my_rank);

    #pragma omp parallel for num_threads(2)
    for (int i = my_rank * Count / size; i < (my_rank + 1) * Count / size; i++)
    {
        for (int k = 0; k < Count; k++)
        {
            if (i != k)
            {                    
                
                //Calculate Distance
                double distance = calculate_distance(stars, i, k);   //Gravitational Distance between 2 Stars

                //Calculate force between two stars - //Gravitational Force between 2 Stars
                double temp_force = (gravitational_constat * (*(stars + i * 4 + 3)) * (*(stars + k * 4 + 3)) / (pow(distance, 2)));
                
                // Add to total
                gravitational_force_list[i - my_rank * Count / size] += temp_force;

                temp_force = 0;
            }
        }
    }

    printf("Calculating was finished by %d.Process!\n", my_rank);
}

//Calculate distance between two stars
double calculate_distance(int *stars, int main_star_index, int other_star_index)
{

    //Get 3 locations and and subtract
    double main_distance = 0;
    double distance_1D = *(stars + main_star_index * 4 + 0) - *(stars + other_star_index * 4 + 0);
    double distance_2D = *(stars + main_star_index * 4 + 1) - *(stars + other_star_index * 4 + 1);
    double distance_3D = *(stars + main_star_index * 4 + 2) - *(stars + other_star_index * 4 + 2);

    //get exponential
    double squared_distances = pow(distance_1D, 2) + pow(distance_2D, 2) + pow(distance_3D, 2);

    //get squared
    main_distance = sqrt(squared_distances);
    return main_distance;
}

//Calculate 2 norm
double calculatePartial_2Norm(int Count, double gravitational_force_list[], int size)
{

    double total = 0;
    int i = 0;

    #pragma omp parallel for num_threads(2) reduction(+:total)
    for (i = 0; i < Count / size; i++)
    {
        total += pow(gravitational_force_list[i], 2);
    }

    return total;
}

//Write to file
void write_file(double gravitational_force_list[], int Count)
{

    FILE *file;
    file = fopen("forces.txt", "w");

    for (int i = 0; i < Count; i++)
    {
        fprintf(file, "%lf \n", gravitational_force_list[i]);
    }

    printf("Gravitational Forces were writed to file! \n");

    fclose(file);
}