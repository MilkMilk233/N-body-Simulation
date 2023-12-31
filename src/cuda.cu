#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int block_size = 512;


int n_body;
int n_iteration;


__global__ void update_position(double *x, double *y, double *vx, double *vy, double *ax, double *ay, int n) {
    //TODO: update position 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
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

__global__ void update_velocity(double *m, double *x, double *y, double *vx, double *vy, double *ax, double *ay, int n) {
    //TODO: calculate force and acceleration, update velocity
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    // printf("Block ID = %d\n",blockIdx.x);
    if (i < n) {  
        double distance, F, v_si, v_si_new, v_sj, v_sj_new;
        for(int j = i+1; j < n; j++){
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



void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];
    double* ax = new double[n_body];
    double* ay = new double[n_body];


    generate_data(m, x, y, vx, vy, ax, ay, n_body);

    Logger l = Logger("cuda", n_body, bound_x, bound_y);

    double *device_m;
    double *device_x;
    double *device_y;
    double *device_vx;
    double *device_vy;
    double *device_ax;
    double *device_ay;

    cudaMalloc(&device_m, n_body*sizeof(double));
    cudaMalloc(&device_x, n_body*sizeof(double));
    cudaMalloc(&device_y, n_body*sizeof(double));
    cudaMalloc(&device_vx, n_body*sizeof(double));
    cudaMalloc(&device_vy, n_body*sizeof(double));
    cudaMalloc(&device_ax, n_body*sizeof(double));
    cudaMalloc(&device_ay, n_body*sizeof(double));

    cudaMemcpy(device_m, m, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_x, x, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, y, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vx, vx, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vy, vy, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_ax, vx, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_ay, vy, n_body*sizeof(double), cudaMemcpyHostToDevice);

    int n_block = n_body / block_size + 1;

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        update_velocity<<<n_block, block_size>>>(device_m, device_x, device_y, device_vx, device_vy, device_ax, device_ay, n_body);
        update_position<<<n_block, block_size>>>(device_x, device_y, device_vx, device_vy, device_ax, device_ay, n_body);

        cudaMemcpy(x, device_x, n_body*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(y, device_y, n_body*sizeof(double), cudaMemcpyDeviceToHost);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        
        // printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

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

    cudaFree(device_m);
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_vx);
    cudaFree(device_vy);
    cudaFree(device_ax);
    cudaFree(device_ay);

    delete m;
    delete x;
    delete y;
    delete vx;
    delete vy;
    delete ax;
    delete ay;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation CUDA Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif

    master();

    printf("Student ID: 120090222\n"); // replace it with your student id
    printf("Name: Chen Zhixin\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation CUDA Implementation\n");

    return 0;

}


