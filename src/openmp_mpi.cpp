/*  
    CSC4005 Project 3 @ Chen Zhixin
    This is the bonus part, combing openMP with MPI.
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;


int my_rank;
int world_size;
int n_omp_threads;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


void update_position(double *x, double *y, double *local_x, double *local_y,
                     double *vx, double *vy, double *ax, double *ay, int start, int end, int i) {
    //TODO: update position 
    vx[i] += ax[i] * dt;
    vy[i] += ay[i] * dt;
    local_x[i] = x[i+start] + vx[i] * dt;
    local_y[i] = y[i+start] + vy[i] * dt;
    // Boundary collision
    while(true){
        if(local_x[i] < sqrt(radius2)){
            local_x[i] += 2*(sqrt(radius2)-local_x[i]);
            vx[i] = -vx[i];
        }
        if(local_x[i] > bound_x-sqrt(radius2)){
            local_x[i] -= 2 * (local_x[i] - (bound_x-sqrt(radius2))); 
            vx[i] = -vx[i];
        }
        if(local_x[i] > 0){
            break;
        }
    }
    while(true){
        if(local_y[i] < sqrt(radius2)){
            local_y[i] += 2*(sqrt(radius2)-local_y[i]);
            vy[i] = -vy[i];
        }
        if(local_y[i] > bound_y-sqrt(radius2)){
            local_y[i] -= 2 * (local_y[i] - (bound_y-sqrt(radius2)));
            vy[i] = -vy[i];
        }
        if(local_y[i] > 0){
            break;
        }
    }
}


void update_velocity(double *m, double *x, double *y, double *vx, double *vy,
                    double *total_vx, double *total_vy, double *ax, double *ay, int start, int end, int i) {
    //TODO: calculate force and acceleration, update velocity
    double distance, F, v_si, v_si_new, v_sj;
    for(int j = 0; j < n_body; j++){
        if(i + start != j){
            distance = sqrt(pow(x[j] - x[i + start], 2) + pow(y[j] - y[i + start], 2));
            if(pow(distance,2) <= 4 * radius2){
                // Elastic collision between two balls
                v_si = (vx[i] * (x[j]-x[i + start]) + vy[i] * (y[j] - y[i + start])) / (distance + err);
                v_sj = (total_vx[j] * (x[i + start]-x[j]) + total_vy[j] * (y[i + start] - y[j])) / (distance + err);
                v_si_new = (v_si * (m[i + start] - m[j]) + 2 * m[j] * v_sj) / (m[i + start] + m[j] + err);
                vx[i] = (v_si_new * (x[j] - x[i + start])) / (distance + err);
                vy[i] = (v_si_new * (y[j] - y[i + start])) / (distance + err);
            }
            else{
                F = ((gravity_const * m[i + start] * m[j]) / (pow(distance,2) + err));
                ax[i] += (F * (x[j] - x[i + start])) / (m[i + start] * distance + err);
                ay[i] += (F * (y[j] - y[i + start])) / (m[i + start] * distance + err);
            }
        }
    }    
}


void slave(){
    MPI_Barrier(MPI_COMM_WORLD);
    // TODO: MPI routine
    int local_size = n_body / world_size;
    int extra_length =  n_body - local_size * world_size;
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];
    double* local_x = new double[local_size+extra_length];
    double* local_y = new double[local_size+extra_length];
    double* local_vx = new double[local_size+extra_length];
    double* local_vy = new double[local_size+extra_length];
    double* local_ax = new double[local_size+extra_length];
    double* local_ay = new double[local_size+extra_length];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(total_m, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < n_iteration; i++){
        // Info retrival
        MPI_Bcast(total_x, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_y, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_vx, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_vy, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(total_vx, local_size, MPI_DOUBLE, local_vx, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(total_vy, local_size, MPI_DOUBLE, local_vy, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(my_rank == world_size - 1 && extra_length){
            MPI_Recv(&local_vx[local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&local_vy[local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Calculation in parallel with openMP
        if(my_rank == world_size-1){
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int j = 0; j < local_size+extra_length; j++) {
                update_velocity(total_m, total_x, total_y, local_vx, local_vy, total_vx, total_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size + extra_length, j);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int j = 0; j < local_size+extra_length; j++) {
                update_position(total_x, total_y, local_x, local_y, local_vx, local_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size + extra_length, j);
            }     
        }
        else{
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int j = 0; j < local_size; j++) {
                update_velocity(total_m, total_x, total_y, local_vx, local_vy, total_vx, total_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size, j);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int j = 0; j < local_size; j++) {
                update_position(total_x, total_y, local_x, local_y, local_vx, local_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size, j);
            }  
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Info gathering
        MPI_Gather(local_vx, local_size, MPI_DOUBLE, total_vx, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_vy, local_size, MPI_DOUBLE, total_vy, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_x, local_size, MPI_DOUBLE, total_x, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_y, local_size, MPI_DOUBLE, total_y, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(my_rank == world_size - 1 && extra_length){
            MPI_Send(&local_vx[local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&local_vy[local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&total_x[(my_rank+1) * local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&total_y[(my_rank+1) * local_size], extra_length , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    delete total_m;
    delete total_x;
    delete total_y;
    delete total_vx;
    delete total_vy;
    delete local_x;
    delete local_y;
    delete local_vx;
    delete local_vy;
    delete local_ax;
    delete local_ay;

}



void master() {
    MPI_Barrier(MPI_COMM_WORLD);
    int local_size = n_body / world_size;
    int extra_length =  n_body - local_size * world_size;
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];
    double* local_x = new double[local_size+extra_length];
    double* local_y = new double[local_size+extra_length];
    double* local_vx = new double[local_size+extra_length];
    double* local_vy = new double[local_size+extra_length];
    double* local_ax = new double[local_size+extra_length];
    double* local_ay = new double[local_size+extra_length];

    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(total_m, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    Logger l = Logger("openmp_with_mpi", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        // Info retrival
        MPI_Bcast(total_x, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_y, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_vx, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(total_vy, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(total_vx, local_size, MPI_DOUBLE, local_vx, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(total_vy, local_size, MPI_DOUBLE, local_vy, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(extra_length){
            MPI_Send(&total_vx[local_size*world_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD);
            MPI_Send(&total_vy[local_size*world_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Calculation in parallel with openMP
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < local_size; i++) {
            update_velocity(total_m, total_x, total_y, local_vx, local_vy, total_vx, total_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size, i);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < local_size; i++) {
            update_position(total_x, total_y, local_x, local_y, local_vx, local_vy, local_ax, local_ay, my_rank * local_size, (my_rank+1) * local_size, i);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Info gathering
        MPI_Gather(local_vx, local_size, MPI_DOUBLE, total_vx, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_vy, local_size, MPI_DOUBLE, total_vy, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_x, local_size, MPI_DOUBLE, total_x, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(local_y, local_size, MPI_DOUBLE, total_y, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(extra_length){
            MPI_Recv(&total_vx[my_rank*local_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&total_vy[my_rank*local_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&total_x[my_rank*local_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&total_y[my_rank*local_size], extra_length , MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        // printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(total_x, total_y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = total_x[i];
            yi = total_y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete total_m;
    delete total_x;
    delete total_y;
    delete local_x;
    delete local_y;
    delete local_vx;
    delete local_vy;
    delete total_vx;
    delete total_vy;
    delete local_ax;
    delete local_ay;
}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI + OpenMP Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 120090222\n"); // replace it with your student id
		printf("Name: Chen Zhixin\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI + OpenMP Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

