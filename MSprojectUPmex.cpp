#include "mex.h" //--This one is required
#include "matrix.h"
#include "math.h"
#include <omp.h>

float sumarray(float *myarray,int length);
inline float minimumvalueoftwo(float a, float b);
inline float maximumvalueoftwo(float a, float b);
float maximumofarray(float *myarray,int length);
inline float myabs(float myvalue);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
// [phi1,phi2] = MSprojectUPmex(phi1,phi2,a,epsilon,hz,max_iter,tol);
    // Input phi1
    double* phi1 = mxGetPr(prhs[0]); 
    // Input phi2
    double* phi2 = mxGetPr(prhs[1]);
    // Number of dimension of phi1
    int numdimphi1 = mxGetNumberOfDimensions(prhs[0]);
    // Dimensions of phi1
    const long unsigned int* sizephi1 = mxGetDimensions(prhs[0]);
    // Input a
    float a = (float)mxGetScalar(prhs[2]);
    // Input epsilon
	float epsilon = (float)mxGetScalar(prhs[3]);
    // Input hz
	float hz = (float)mxGetScalar(prhs[4]);
    // Input max_iter
	int max_iter = (int)mxGetScalar(prhs[5]);
    // Input tol
	float tol = (float)mxGetScalar(prhs[6]);
    // Output phi1
    plhs[0] = mxCreateNumericArray(numdimphi1, sizephi1, mxDOUBLE_CLASS, mxREAL);
    double* phi1outReturn = mxGetPr(plhs[0]);     
    // Output phi2
    plhs[1] = mxCreateNumericArray(numdimphi1, sizephi1, mxDOUBLE_CLASS, mxREAL);
    double* phi2outReturn = mxGetPr(plhs[1]);
    //--------------------------------------------------------------------
    // Function
    int lengthTmp1 = (sizephi1[2]+1)*sizephi1[2]*(sizephi1[2]+1);
    int lengthTmp2 = sizephi1[0]*sizephi1[1]*sizephi1[2];
    float* phi1out = new float[sizephi1[0]*sizephi1[1]*(sizephi1[2]+1)];
    float* phi2out = new float[sizephi1[0]*sizephi1[1]*(sizephi1[2]+1)];
    // Initialize output
	for(int i=0; i<lengthTmp2; ++i) {
		phi1out[i] = (float)phi1[i];
		phi2out[i] = (float)phi2[i];
	}
    // Initialize variables
	const long unsigned int sizelastDiffPhi1[3] = {sizephi1[2]+1,sizephi1[2],sizephi1[2]+1};
    float* lastDiffPhi1 = new float[lengthTmp1];
    float* lastDiffPhi2 = new float[lengthTmp1];
    float* lastDiffPhi1Add = new float[sizephi1[2]];
    float* lastDiffPhi2Add = new float[sizephi1[2]];
	float* aux1 = new float[sizephi1[2]];
	float* aux2 = new float[sizephi1[2]];
	float* phi1old = new float[sizephi1[2]];
	float* phi2old = new float[sizephi1[2]];
	float* maxarray1 = new float[sizephi1[2]];
	float* maxarray2 = new float[sizephi1[2]];
    // Precompute DeltaS and powers
	float* powers = new float[sizephi1[2]+1];
	float* DeltaS = new float[sizephi1[2]+1];
    for(int i=0; i < sizephi1[2]+1; i++) {
        DeltaS[i] = hz * (float)i;
    }
    int indexTmp1 = 0; int indexTmp2 = 0;
    float A = 0; float B = 0; float beta = 0;
    float norm = 0; float C = 0; float factor;
    int iter; float diff;
    // Iteration for every pixel
    for(int pixelx=0; pixelx<sizephi1[0]; pixelx++) {
        for(int pixely=0; pixely<sizephi1[1]; pixely++) {
			int lengthTmp = lengthTmp1;
            // Set lastDiff to zero
            for(int i=0; i < lengthTmp; i++) {
				lastDiffPhi1[i] = 0.0f;
				lastDiffPhi2[i] = 0.0f;
			}
            // Set iteration counter and error 
            int iter = 0; float error = 1.0f;
            // Iteration
            while( iter <= max_iter && error > tol ) {   
                error = 0.0f;
                for(int t1=0; t1<sizephi1[2]; ++t1) {
                    for(int t2=t1; t2<sizephi1[2]; ++t2) {                         
						lengthTmp = t2-t1+1;
                        // Undo last projection
                        for(int i=0; i < lengthTmp; ++i) {
                            indexTmp1 = pixelx + pixely*sizephi1[0] + (t1+i)*sizephi1[0]*sizephi1[1];
                            aux1[i] = phi1out[indexTmp1];
                            aux2[i] = phi2out[indexTmp1];
                            indexTmp2 = t1 + i*sizelastDiffPhi1[0] + t2*sizelastDiffPhi1[0]*sizelastDiffPhi1[1];
                            phi1old[i] = aux1[i] - lastDiffPhi1[indexTmp2];
                            phi2old[i] = aux2[i] - lastDiffPhi2[indexTmp2];
                        }
                        // Compute (A,B) and beta
                        A = hz * sumarray(phi1old,lengthTmp); // Integral phi1
                        B = hz * sumarray(phi2old,lengthTmp); // Integral phi2
                        norm = sqrt(A*A+B*B)+1e-7;
                        C = minimumvalueoftwo(DeltaS[lengthTmp]+epsilon,a*DeltaS[lengthTmp]);
                        beta = minimumvalueoftwo(0.0f, C/(norm*DeltaS[lengthTmp])-1.0f/DeltaS[lengthTmp]);
                        // Project phiold
                        for(int i=0; i<lengthTmp; ++i) {
                            indexTmp1 = pixelx + pixely*sizephi1[0] + (t1+i)*sizephi1[0]*sizephi1[1];
                            phi1out[indexTmp1] = phi1old[i] + beta * A;
                            phi2out[indexTmp1] = phi2old[i] + beta * B;
                        }
                        // Update lastDiffPhi
                        for(int s=0; s<lengthTmp; ++s) {
                            indexTmp1 = t1 + s*sizelastDiffPhi1[0] + t2*sizelastDiffPhi1[0]*sizelastDiffPhi1[1];
                            indexTmp2 = pixelx + pixely*sizephi1[0] + (t1+s)*sizephi1[0]*sizephi1[1];
                            lastDiffPhi1[indexTmp1] = phi1out[indexTmp2] - phi1old[s];
                            lastDiffPhi2[indexTmp1] = phi2out[indexTmp2] - phi2old[s];
                        }
                        // Update error
                         for(int i=0; i<lengthTmp; ++i) {
                             indexTmp1 = pixelx + pixely*sizephi1[0] + (t1+i)*sizephi1[0]*sizephi1[1];
                             maxarray1[i] = myabs(aux1[i] - phi1out[indexTmp1]);
                             maxarray2[i] = myabs(aux2[i] - phi2out[indexTmp1]);
                         }
                         error = maximumvalueoftwo(diff, maximumvalueoftwo(maximumofarray(maxarray1,lengthTmp), maximumofarray(maxarray2,lengthTmp)));
                    }
                }
                //----------------------------------------------------
                // Projection onto additional constraint |phi|<=a
                for(int t1=0; t1<sizephi1[2]; ++t1) {
                    indexTmp1 = sizephi1[2] + t1*sizelastDiffPhi1[0] + sizephi1[2]*sizelastDiffPhi1[0]*sizelastDiffPhi1[1];
                    indexTmp2 = pixelx + pixely*sizephi1[0] + t1*sizephi1[0]*sizephi1[1];
                    // Undo last projection
                    phi1old[t1] = phi1out[indexTmp2] - lastDiffPhi1[indexTmp1];
                    phi2old[t1] = phi2out[indexTmp2] - lastDiffPhi2[indexTmp1];
                    // Projection
                    norm = sqrt(phi1old[t1]*phi1old[t1]+phi2old[t1]*phi2old[t1])+1e-7;
                    factor = minimumvalueoftwo ( 1.0f, a/norm );
                    phi1out[indexTmp2] = phi1old[t1]*factor;
                    phi2out[indexTmp2] = phi2old[t1]*factor;
                    // Save lastDiffPhi
                    lastDiffPhi1[indexTmp1] = phi1out[indexTmp2] - phi1old[t1];
                    lastDiffPhi2[indexTmp1] = phi2out[indexTmp2] - phi2old[t1];
                }    
                //----------------------------------------------------
				iter++;
            }
        }
    }
    // Return values 
    for(int i=0; i<lengthTmp2; ++i) {
		phi1outReturn[i] = (double)phi1out[i];
		phi2outReturn[i] = (double)phi2out[i];
	}
    // Delete variables
    delete [] phi1out;
    delete [] phi2out;
    delete [] lastDiffPhi1;
    delete [] lastDiffPhi2;
    delete [] aux1;
    delete [] aux2;
    delete [] phi1old;
    delete [] phi2old;
    delete [] maxarray1;
    delete [] maxarray2;
    delete [] powers;
    delete [] DeltaS;
    return;
}

// Sum of all entries of an array
float sumarray(float *myarray, int length) {
    float mysum = 0;
    for(int i=0; i<length; i++) mysum += myarray[i];
    return mysum;
}

// Minimum of two doubles
inline float minimumvalueoftwo(float a, float b) {
    return a > b ? b : a;
}

// Maximum of two doubles
inline float maximumvalueoftwo(float a, float b) {
   return a < b ? b : a;
}

// Maximum entry of an array
float maximumofarray(float *myarray, int length) {
    float maximum = myarray[0];
    for(int i=1; i<length; i++) {
        if (myarray[i] > maximum) maximum = myarray[i];
    }
    return maximum;
}

// Absolute value
inline float myabs(float myvalue) {
    return myvalue > 0 ? myvalue : -myvalue;
}
