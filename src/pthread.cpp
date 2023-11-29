#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;
int ready_thread;
pthread_mutex_t count_mutex;
pthread_cond_t count_threshold_cv;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, double *ax, double *ay, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
        ax[i] = 0.0f;
        ay[i] = 0.0f;
    }
}

void update_position(double *x, double *y, double *vx, double *vy, double *ax, double *ay, int start, int end) {
    
    // Equivalent to MPI_Barrier in MPI version, in case of data race
    pthread_mutex_lock(&count_mutex);
    ready_thread += 1;
    if(ready_thread < n_thd){
        pthread_cond_wait(&count_threshold_cv, &count_mutex);
    }
    else{
        pthread_cond_broadcast(&count_threshold_cv);
    }
    ready_thread = 0;
    pthread_mutex_unlock(&count_mutex);

    //TODO: update position 
    for(int i = start; i < end; i++){
        vx[i] += ax[i] * dt;
        vy[i] += ay[i] * dt;
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        // Boundary collision
        while(true){
            if(x[i] < sqrt(radius2)){
                x[i] += 2*(sqrt(radius2)-x[i]);
                vx[i] = -vx[i];
            }
            if(x[i] > bound_x-sqrt(radius2)){
                x[i] -= 2 * (x[i] - (bound_x-sqrt(radius2))); 
                vx[i] = -vx[i];
            }
            if(x[i] > 0){
                break;
            }
        }
        while(true){
            if(y[i] < sqrt(radius2)){
                y[i] += 2*(sqrt(radius2)-y[i]);
                vy[i] = -vy[i];
            }
            if(y[i] > bound_y-sqrt(radius2)){
                y[i] -= 2 * (y[i] - (bound_y-sqrt(radius2)));
                vy[i] = -vy[i];
            }
            if(y[i] > 0){
                break;
            }
        }
    }

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, double *ax, double *ay, int start, int end) {
    //TODO: calculate force and acceleration, update velocity
    double distance, F, v_si, v_si_new, v_sj, v_sj_new;
    for(int i = start; i < end ; i++){
        for(int j = i+1; j < n_body; j++){
            distance = sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2));
            if(pow(distance,2) <= 4 * radius2){
                // Elastic collision between two balls
                v_si = (vx[i] * (x[j]-x[i]) + vy[i] * (y[j] - y[i])) / (distance + err);
                v_sj = (vx[j] * (x[i]-x[j]) + vy[j] * (y[i] - y[j])) / (distance + err);
                v_si_new = (v_si * (m[i] - m[j]) + 2 * m[j] * v_sj) / (m[i] + m[j] + err);
                v_sj_new = (v_sj * (m[j] - m[i]) + 2 * m[i] * v_si) / (m[i] + m[j] + err);
                vx[i] = (v_si_new * (x[j] - x[i])) / (distance + err);
                vy[i] = (v_si_new * (y[j] - y[i])) / (distance + err);
                vx[j] = (v_sj_new * (x[i] - x[j])) / (distance + err);
                vy[j] = (v_sj_new * (y[i] - y[j])) / (distance + err);
            }
            else{
                F = ((gravity_const * m[i] * m[j]) / (pow(distance,2) + err));
                ax[i] += (F * (x[j] - x[i])) / (m[i] * distance + err);
                ay[i] += (F * (y[j] - y[i])) / (m[i] * distance + err);
                ax[j] += (F * (x[i] - x[j])) / (m[j] * distance + err);
                ay[j] += (F * (y[i] - y[j])) / (m[j] * distance + err);
            }
        }
    }
}


typedef struct {
    //TODO: specify your arguments for threads
    int start;
    int end;
    double* m;
    double* x;
    double* y;
    double* vx;
    double* vy;
    double* ax;
    double* ay;
    //TODO END
} Args;


void* worker(void* args) {
    //TODO: procedure in each threads
    Args* my_arg = (Args*) args;
    int start = my_arg->start;
    int end = my_arg->end;
    double* m = my_arg->m;
    double* x = my_arg->x;
    double* y = my_arg->y;
    double* vx = my_arg->vx;
    double* vy = my_arg->vy;
    double* ax = my_arg->ax;
    double* ay = my_arg->ay;
    update_velocity(m, x, y, vx, vy, ax, ay, start, end);
    update_position(x, y, vx, vy, ax, ay, start, end);
    pthread_exit(NULL);
    // TODO END
}


void master(){
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];
    double* ax = new double[n_body];
    double* ay = new double[n_body];

    generate_data(m, x, y, vx, vy, ax, ay, n_body);

    pthread_mutex_init(&count_mutex, NULL);
    pthread_cond_init(&count_threshold_cv, NULL);
    pthread_t thds[n_thd]; // thread pool
    Args args[n_thd]; // arguments for all threads
    int my_data_size = n_body / n_thd;
    int current_size = 0;
    // Configure arguments for each thread
    for (int thd = 0; thd < n_thd; thd++){
        args[thd].start = current_size;
        args[thd].m = m;
        args[thd].x = x;
        args[thd].y = y;
        args[thd].vx = vx;
        args[thd].vy = vy;
        args[thd].ax = ax;
        args[thd].ay = ay;
        if(thd != n_thd-1){
            current_size += my_data_size;
            args[thd].end = current_size;
        }
        else{
            args[thd].end = n_body;
        }
    }

    ready_thread = 0;   // Counting tool for "barrier"
    Logger l = Logger("pthread", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs
        // Create threads in every single iteration
        for (int thd = 0; thd < n_thd; thd++) pthread_create(&thds[thd], NULL, worker, &args[thd]);
        for (int thd = 0; thd < n_thd; thd++) pthread_join(thds[thd], NULL);
        //TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);
        l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete m;
    delete x;
    delete y;
    delete vx;
    delete vy;
    delete ax;
    delete ay;
    pthread_mutex_destroy(&count_mutex);
	pthread_cond_destroy(&count_threshold_cv);

}


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

	return 0;
}

