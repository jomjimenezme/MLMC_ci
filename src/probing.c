#include "main.h"
#include "data_layout.h"

int max(int v[], int size) {
    int max = v[0];  

    for (int i = 1; i < size; i++) { 
        if (v[i] > max) {          
            max = v[i];
        }
    }

    return max;
}

void vector_copy(int *dest, int *src, int size) {
    for (int i = 0; i < size; i++) {
        dest[i] = src[i];
    }
}

/*
void colors_init(){
    MALLOC(g.num_colors, int, g.num_levels);
    
    g.colors = (int**)malloc(g.num_levels * sizeof(int*));
    
    if (g.colors == NULL)
        error0("Allocation error0\n");
}
*/

void setup_local_colors(){
    
    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    
    g.local_colors = (int**)malloc(g.num_levels * sizeof(int*));
    
    for(int level = 0; level < g.num_levels; level++){
        
        int T = g.global_lattice[level][0];
        int Z = g.global_lattice[level][1];
        int Y = g.global_lattice[level][2];
        int X = g.global_lattice[level][3];
        
        int size = T * Z * Y * X;
        int local_size = size/num_processes;
        int *current_level_colors;
        
        if(g.my_rank == 0){
            MALLOC(current_level_colors, int, size);
            vector_copy(current_level_colors, g.colors[level], size);
        }
        
        MALLOC(g.local_colors[level], int, local_size);
    
        int *current_level_local_colors;
        MALLOC(current_level_local_colors, int, local_size);
    
        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Scatter(current_level_colors, local_size, MPI_INT,
                current_level_local_colors, local_size, MPI_INT, 
                0, MPI_COMM_WORLD);
    
        vector_copy(g.local_colors[level], current_level_local_colors, local_size);
        
        FREE(current_level_colors, int, size);
        FREE(current_level_local_colors, int, local_size);
    }
    
    if(g.my_rank==0){
        for(int i = 0; i < g.num_levels; i++){
	  int T = g.global_lattice[i][0];
	  int Z = g.global_lattice[i][1];
	  int Y = g.global_lattice[i][2];
	  int X = g.global_lattice[i][3];
          int size = T * Z * Y * X;
          FREE(g.colors[i], int*, size );
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(g.num_colors, g.num_levels, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   
}

void coloring_scheme(){
    
//colors_init();

for(int level = 0; level < g.num_levels; level++){
    
    int T = g.global_lattice[level][0];
    int Z = g.global_lattice[level][1];
    int Y = g.global_lattice[level][2];
    int X = g.global_lattice[level][3];
    
    int total_points = T * Z * Y * X;
    
    int size[4];
    
    size[0] = T;
    size[1] = Z;
    size[2] = Y;
    size[3] = X;
    
    int coordinates[4];
    
    MALLOC(g.colors[level], int, total_points);
    
    for(int t = 0; t < T; t++){
        coordinates[0] = t;
        for(int z = 0; z < Z; z++){
            coordinates[1] = z;
            for(int y = 0; y < Y; y++){
                coordinates[2] = y;
                for(int x = 0; x < X; x++){
                    coordinates[3] = x;
                    int idx = lex_index(t, z, y, x, size);

                    int col = 0;
                    int power = 1;
                    for(int k = 0; k<4; k++){
                        int w_tilde = coordinates[k]%(g.coloring_distance + 1);
                        col += w_tilde * power;
                        power *= (g.coloring_distance + 1);
	               }    
	               g.colors[level][idx] = col + 1;
               }
            }
        }
    }
    g.num_colors[level] = max(g.colors[level], total_points);
}

}


void generate_neighbors(int t, int z, int y, int x, int **neighbors, int *num_neighbors, int size[4]) {
    int T = size[0];
    int Z = size[1];
    int Y = size[2];
    int X = size[3];
    
    int max_neighbors = 1; // Includes central point
    for (int delta = 1; delta <= g.coloring_distance; delta++) {
        max_neighbors += 8 * delta * delta * delta; // Estimate maximum number
    }
    *neighbors = (int *)malloc(max_neighbors * sizeof(int));
    *num_neighbors = 0;

    // Iterate over all displacements combination
    for (int dt = -g.coloring_distance; dt <= g.coloring_distance; dt++) {
        for (int dz = -g.coloring_distance; dz <= g.coloring_distance; dz++) {
            for (int dy = -g.coloring_distance; dy <= g.coloring_distance; dy++) {
                for (int dx = -g.coloring_distance; dx <= g.coloring_distance; dx++) {
                    
                    //Verify that the distance is the right one
                    if (abs(dt) + abs(dz) + abs(dy) + abs(dx) > g.coloring_distance) {
                        continue;
                    }

                    // Skip central point (0, 0, 0, 0)
                    if (dt == 0 && dz == 0 && dy == 0 && dx == 0) {
                        continue;
                    }

                    // Compute coordinates of the neighbor
                    int t_neighbor = (t + dt + T) % T;
                    int z_neighbor = (z + dz + Z) % Z;
                    int y_neighbor = (y + dy + Y) % Y;
                    int x_neighbor = (x + dx + X) % X;

                    // Compute the lexicographic index of the neighbor
                    (*neighbors)[*num_neighbors] = lex_index(t_neighbor, z_neighbor, y_neighbor, x_neighbor, size);
                    (*num_neighbors)++;
                }
            }
        }
    }
}


void graph_coloring() {
    
    MALLOC(g.num_colors, int, g.num_levels);
    
    if(g.my_rank == 0){
        
    double time_taken;
    
    double start_time = MPI_Wtime();
    
    g.colors = (int**)malloc(g.num_levels * sizeof(int*));
    
    if (g.colors == NULL)
        error0("Allocation error0\n");
    
    for(int level = 0; level < g.num_levels; level++){
    
    int T = g.global_lattice[level][0];
    int Z = g.global_lattice[level][1];
    int Y = g.global_lattice[level][2];
    int X = g.global_lattice[level][3];
    
    int size[4];
    
    size[0] = T;
    size[1] = Z;
    size[2] = Y;
    size[3] = X;

    int total_points = T * Z * Y * X;
    
    MALLOC(g.colors[level], int, total_points);

    // Set all colors to 0 (not assigned)
    for (int i = 0; i < total_points; i++) {
        g.colors[level][i] = 0;
    }

    // Iterate over all lattice sites
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < Z; z++) {
            for (int y = 0; y < Y; y++) {
                for (int x = 0; x < X; x++) {
                    int index = lex_index(t, z, y, x, size);

                    // Skip if the site has already been assigned a color
                    if (g.colors[level][index] != 0) {
                        continue;
                    }

                    // Find neighbors within distance d
                    int *neighbors;
                    int num_neighbors;
                    generate_neighbors(t, z, y, x, &neighbors, &num_neighbors, size);

                    int used_colors[256] = {0}; 
                    for (int i = 0; i < num_neighbors; i++) {
                        int neighbor_index = neighbors[i];
                        if (g.colors[level][neighbor_index] != 0) {
                            used_colors[g.colors[level][neighbor_index] - 1] = 1;
                        }
                    }

                    // Assign the minimum possible color
                    int color = 1;
                    while (used_colors[color - 1]) {
                        color++;
                    }
                    g.colors[level][index] = color;

                    // Free memory of neighbors
                    free(neighbors);
                }
            }
        }
    }
    g.num_colors[level] = max(g.colors[level], total_points);
    
    }
    
    double end_time = MPI_Wtime();
    
    time_taken = end_time - start_time;
                 
    printf("\nTime for coloring: %f seconds\n", time_taken);
    
    /*
    FILE *file = fopen("/home/papace/c_code/mpi_parallelization/print_files/colors.txt", "w");
    if (file == NULL) {
        perror("Errore nell'apertura del file");
    }
    
    for(int i = 0; i < g.num_levels; i++){
        int T = g.global_lattice[i][0];
    int Z = g.global_lattice[i][1];
    int Y = g.global_lattice[i][2];
    int X = g.global_lattice[i][3];
        int size = T * Z * Y * X;
        printf("size = %d",size);
        fprintf(file, "\ncolors at level %d\n[", i);
        for(int j = 0; j < size; j ++){
            fprintf(file, " %d ", g.colors[i][j]);
        }
        fprintf(file, "]\n");
    }
    fclose(file);
    */
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    setup_local_colors();
    //setup_num_colors();
}


// DOWN FOR A MIXTURE OF THE COLORING METHODS 
// UP FOR JUST ONE OF THE TWO ---> generate_neighbors() function is needed for both, do not comment it!

/*
void greedy_coloring(int *colors, int total_points, int size[4]) {

    // Set all colors to 0 (not assigned)
    for (int i = 0; i < total_points; i++) {
        colors[i] = 0;
    }

    // Iterate over all lattice sites
    for (int t = 0; t < size[0]; t++) {
        for (int z = 0; z < size[1]; z++) {
            for (int y = 0; y < size[2]; y++) {
                for (int x = 0; x < size[3]; x++) {
                    int index = lex_index(t, z, y, x, size);

                    // Skip if the site has already been assigned a color
                    if (colors[index] != 0) {
                        continue;
                    }

                    // Find neighbors within distance d
                    int *neighbors;
                    int num_neighbors;
                    generate_neighbors(t, z, y, x, &neighbors, &num_neighbors, size);

                    int used_colors[256] = {0}; 
                    for (int i = 0; i < num_neighbors; i++) {
                        int neighbor_index = neighbors[i];
                        if (colors[neighbor_index] != 0) {
                            used_colors[colors[neighbor_index] - 1] = 1;
                        }
                    }

                    // Assign the minimum possible color
                    int color = 1;
                    while (used_colors[color - 1]) {
                        color++;
                    }
                    colors[index] = color;

                    // Free memory of neighbors
                    free(neighbors);
                }
            }
        }
    }
    
    }

void coloring_scheme(int *colors, int size[4]){
    
    int coordinates[4];
    
    for(int t = 0; t < size[0]; t++){
        coordinates[0] = t;
        for(int z = 0; z < size[1]; z++){
            coordinates[1] = z;
            for(int y = 0; y < size[2]; y++){
                coordinates[2] = y;
                for(int x = 0; x < size[3]; x++){
                    coordinates[3] = x;
                    int idx = lex_index(t, z, y, x, size);

                    int col = 0;
                    int power = 1;
                    for(int k = 0; k<4; k++){
                        int w_tilde = coordinates[k]%(g.coloring_distance + 1);
                        col += w_tilde * power;
                        power *= (g.coloring_distance + 1);
	               }    
	               colors[idx] = col + 1;
               }
            }
        }
    }
    
}



void graph_coloring(){
    
    struct timespec start, end;
    double time_taken;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    colors_init();
    
    for(int level = 0; level < g.num_levels; level++){
        
        int T = g.global_lattice[level][0];
        int Z = g.global_lattice[level][1];
        int Y = g.global_lattice[level][2];
        int X = g.global_lattice[level][3];
    
        int size[4];
    
        size[0] = T;
        size[1] = Z;
        size[2] = Y;
        size[3] = X;

        int total_points = T * Z * Y * X;
    
        MALLOC(g.colors[level], int, total_points);
        
        if(level < g.coloring_method)
            coloring_scheme(g.colors[level], size);
        else
            greedy_coloring(g.colors[level], total_points, size);
        
        g.num_colors[level] = max(g.colors[level], total_points);
        
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    time_taken = (end.tv_sec - start.tv_sec) +
                 (end.tv_nsec - start.tv_nsec) / 1e9;
                 
    if(g.my_rank == 0)
        printf("\nTime for coloring: %f seconds\n", time_taken);
    
}
*/
