
/* This is our assignment for the SIR model */


/* Including various libraries */
#include <iostream>
#include <vector>
#include <fstream>

/* Defining constants */
const float N=1000;
const float beta = 0.2;
const float gamma = 1.0/10.0;

/* Defining or constant step size, without adaptive dt because we are n00bs */
const float t_step = 0.01;

/* Defining the number of iterations, i.e. days that we will run the model for */
const int N_ite = 1e5;

/* Initial parameters */
std::vector<float> I = {1.0};
std::vector<float> S = {N-I.back()};
std::vector<float> R = {0.0};

/* Define a function calculating next step of S_(n+1) */
void update_vectors( std::vector<float> &S, std::vector<float> &I, std::vector<float> &R  ) {
    
    /* Take the previous S value, aka the last element in the vector. Same for I and R */
    float S_now = S.back();
    float I_now = I.back();
    float R_now = R.back();
    
    /* Discretized equation for dS ------------------------ */
    float S_next = S_now - ( beta * I_now * (S_now/N) )*t_step;
        
    /* Append S_next to the S vector */
    S.push_back(S_next);
    
    
    /* Discretized equation for dI ------------------------ */
    float I_next = I_now + ( beta * I_now * (S_now/N) - gamma * I_now )*t_step;
    
    /* Append I_next to the I vector */
    I.push_back(I_next);
    
    
    /* Discretized equation for dR ------------------------ */
    float R_next = R_now + ( gamma * I_now )*t_step;
    
    /* Append I_next to the I vector */
    R.push_back(R_next);  
}

/* Update the vector N_ite times, so we will have 1+N_ite long vectors with all our values */
int main() {
    
    /* Define name of file with results */
    std::ofstream myfile;
    std::string fname = "./sir_output.txt";
    
    /* Open the file so we can write to it */
    myfile.open(fname);
    
    /* Define the header */
    myfile << "SIR \n";
    
    /* Insert the initial parameters on the first row */
    myfile << S.back() << ' '
           << I.back() << ' '
           << R.back() << '\n';
    
    /* Make a for loop to update the vectors and export values to file ----------- */
    for (int i = 0; i < N_ite; ++i) {  
        update_vectors(S, I, R);
        
        /* Insert to file, a row for each time step containing S, I, R */
        myfile << S.back() << ' '
               << I.back() << ' '
               << R.back() << '\n';
    }

    /* Close the file */
    myfile.close();
    
    return 0;
}

