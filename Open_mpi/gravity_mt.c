#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

void Read_file(char* File_name, int Count ,int *stars);
void calculate_gravitational_force(int Count ,int *stars,double gravitational_force_list[],int treat_count);
double calculate_distance(int *stars,int main_star_index,int other_star_index);
void calculate_2Norm(int Count,double gravitational_force_list[],int thread_count);
void write_file(double gravitational_force_list[],int Count);


int main(int argc, char* argv[]){

    char* File_name;
    int Star_count;
    int thread_count;

    // get variables from input
    File_name=argv[1];
    Star_count=atoi(argv[2]);
    thread_count=atoi(argv[3]);
    
    //Allocate stars array pointer and graviatational force array with malloc
    int *stars = (int *)malloc(Star_count * 4 * sizeof(int));                           //Star İnformation Array
    double *gravitational_force_list = (double *)malloc(Star_count * sizeof(double));   //Gravitational Force List
    
    //Get Stars From File
    Read_file(File_name,Star_count,stars);

    //Calculate Gravitational Force
    calculate_gravitational_force(Star_count ,stars,gravitational_force_list,thread_count);

    //Calculate 2 Norm
    calculate_2Norm(Star_count,gravitational_force_list,thread_count);

    //Write to file
    write_file(gravitational_force_list,Star_count);

    free(stars);
    free(gravitational_force_list);

    return 0;
}

//Get Stars From File
void Read_file(char* File_name, int Count ,int *stars){

    FILE *file;
    file = fopen(File_name, "r");

    if(file){

        for(int i=0;i<Count;i++){
            for(int j=0;j<4;j++){
            fscanf(file, "%d", (stars +i*4+j));
            }
        }
        printf("Star informations were got from file ! \n");
    }
    else{
        printf("File can not be opened !! \n");
    }
    fclose(file);
}

//Calculate Fravitational force
void calculate_gravitational_force(int Count ,int *stars,double gravitational_force_list[],int thread_count){

    const double gravitational_constat=0.0006674;    //Gravitational Constant (G)
    double temp_force=0;                             //Gravitational Force between 2 Stars
    double distance=0;                               //Gravitational Distance between 2 Stars
    double start;                                    //Start Calculating time
    double end;                                      //Finish Calculating time                       

    //Filling the pointer inside
    for(int i=0; i<Count; i++){
        gravitational_force_list[i]=0;
    }
    
    printf("Calculating...  ( the process may take a long time ) \n");
    //Start time
    start =  omp_get_wtime();

    # pragma omp parallel for num_threads(thread_count) 
    for(int i=0; i<Count; i++){
        for(int k=0;k<Count;k++){
            if(i!=k){
                //Calculate Distance
                distance= calculate_distance(stars,i,k);
                //Calculate force between two stars
                temp_force=(gravitational_constat * (*(stars +i*4+3)) * (*(stars +k*4+3)) / (pow(distance,2)));
                // Add to total
                gravitational_force_list[i]+=temp_force;
                temp_force=0;
            } 
        }
    }
    // Finish time
    end = omp_get_wtime();
    double time_spent = (double)(end - start);
    printf("Calculating was finished!\n");
    printf("Calculation Time : %lf Seconds \n",time_spent);

}

//Calculate distance between two stars
double calculate_distance(int *stars,int main_star_index,int other_star_index){

    //Get 3 locations and and subtract
    double main_distance=0;
    double distance_1D=*(stars +main_star_index*4+0)- *(stars +other_star_index*4+0);
    double distance_2D=*(stars +main_star_index*4+1)- *(stars +other_star_index*4+1);
    double distance_3D=*(stars +main_star_index*4+2)- *(stars +other_star_index*4+2);

    //get exponential
    double squared_distances=pow(distance_1D,2) + pow(distance_2D,2) + pow(distance_3D,2);

    //get squared
    main_distance= sqrt(squared_distances);
    return main_distance;
}

//Calculate 2 norm
void calculate_2Norm(int Count,double gravitational_force_list[],int thread_count){

    double total=0;
    int i = 0;
     # pragma omp parallel  for num_threads(thread_count)  reduction(+:total)
    for(i=0;i<Count;i++){
        
        # pragma omp critical 
        total+=pow(gravitational_force_list[i],2);
    }
    total= sqrt(total);

    printf("2-Norm :  %2.13E \n",total);
}

//Write to file
void write_file(double gravitational_force_list[],int Count){

    FILE *file;
    file = fopen("forces.txt", "w");

    for(int i = 0; i < Count;i++){
       fprintf (file, "%lf \n",gravitational_force_list[i]);
    }

    printf("Gravitational Forces were writed to file! \n");

    fclose(file);

}