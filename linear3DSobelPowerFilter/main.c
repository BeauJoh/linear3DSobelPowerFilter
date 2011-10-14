//
//  main.c
//  fastFourierTransform
//
//  Created by Beau Johnston on 19/09/11.
//  Copyright 2011 University Of New England. All rights reserved.
//

#include "global.h"
#include "RGBAUtilities.h"
#include "FileHandler.h"

/* Utility functions not needed on Vanilla Essence C kernel */
int width, height;
uint8* readImage(char* fileName){
    readPngFile(fileName);
    
    width = getImageWidth();
    height = getImageLength();
    
    uint8 * buffer = malloc(sizeof(uint8)*getImageSize());    
    memcpy(buffer, getImage(), getImageSize());
    return buffer;
}

void saveImage(char* fileName, uint8*buffer){
    setImage(buffer);
    writePngFile(fileName);
    return;
}


/* Utility functions essential for computation on the Vanilla Essence C kernel */

/*
 This computes an in-place complex-to-complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir =  1 gives forward transform
 dir = -1 gives reverse transform
 */
int FFT(int dir,long m,float *x,float *y){
	long n,i,i1,j,k,i2,l,l1,l2;
	double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
	// Calculate the number of points
	n = 1;
	for (i=0;i<m;i++)
		n *= 2;
    
	// Do the bit reversal
	i2 = n >> 1;
	j = 0;
	for (i=0;i<n-1;i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = (float)tx;
			y[j] = (float)ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
    
	// Compute the FFT
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<n;i+=l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = (float)(x[i] - t1);
				y[i1] = (float)(y[i] - t2);
				x[i] += (float)t1;
				y[i] += (float)t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == FFT_FORWARD)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}
    
	// Scaling for forward transform
	if (dir == FFT_FORWARD) {
		for (i=0;i<n;i++) {
			x[i] /= n;
			y[i] /= n;
		}
	}
    
	return 1;
}

/*
 Direct fourier transform
 */
int DFT(int dir,int m,float *x1,float *y1)
{
    long i,k;
    float arg;
    float cosarg,sinarg;
    float *x2=NULL,*y2=NULL;
    
    x2 = malloc(m*sizeof(float));
    y2 = malloc(m*sizeof(float));
    if (x2 == NULL || y2 == NULL)
        return(FALSE);
    
    for (i=0;i<m;i++) {
        x2[i] = 0;
        y2[i] = 0;
        arg = - dir * 2.0 * 3.141592654 * (float)i / (float)m;
        for (k=0;k<m;k++) {
            cosarg = cos(k * arg);
            sinarg = sin(k * arg);
            x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
            y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
        }
    }
    
    /* Copy the data back */
    if (dir == 1) {
        for (i=0;i<m;i++) {
            x1[i] = x2[i] / (float)m;
            y1[i] = y2[i] / (float)m;
        }
    } else {
        for (i=0;i<m;i++) {
            x1[i] = x2[i];
            y1[i] = y2[i];
        }
    }
    
    free(x2);
    free(y2);
    return(TRUE);
}

//long ClosestPower2(long x){
//	long double temp=log(x)/log(2);
//	return (long)pow(2,((int)(temp+0.5)));
//}

int ClosestPower2(int x){
	//long double temp=log(x)/log(2);
	return log(x)/log(2);
}

int main (int argc, const char * argv[])
{
#define DoingRealWork   
#ifdef DoingRealWork
    
    //char* fileName="../../../../../dataResources/High-Res-Stage-24-Take-4/out.png";
    //char* outputImageFileName = "../../../../../dataResources/output/result.png";
    
    //a copy in executables neighbour used for debugging!
    char* fileName="dataSet/out.png";
    char* outputImageFileName = "dataSetOut/result.png";
    
    generateListOfAssociatedFiles(fileName);
        
    int depth = numberOfFiles();
    
    uint8*bigBuffer;
    bool firstRun = true;
    printf("Set up ... about to read image stack\n");
    system("pwd");
    printf("\n");
    
    //load all images into a buffer
    for (int i = 0; i < numberOfFiles(); i++) {
        char* tmp = getNextFileName();
        //printf("next file name from main is :%s \n", tmp);
        readPngFile(tmp);
        width = getImageWidth();
        height = getImageLength();
        uint8 *buffer = malloc(sizeof(uint8)*getImageSize());
        memcpy(buffer, getImage(), getImageSize());
        if (firstRun) {
            //if its the first run we don't know the dimensions of the image
            //and thus don't know how much memory to statically allocate
            bigBuffer = malloc(sizeof(uint8)*getImageSize()*depth);
            firstRun = false;
        }
        memcpy(bigBuffer+(i*getImageSize()), buffer, getImageSize());
        free(buffer);
    } 
    
    
    /* ----------------------------------------> KERNEL BEGINS HERE! <---------------------------------------- */
    
    int stackSize = sizeof(float)*getImageSize()*numberOfFiles();
    
    //Create Nx * Ny * Nz array of data. i.e. embryo slices (Nx * Ny) pixels * (Nz) slices deep. Denoted by Da.
    float * DaR = normaliseStack(bigBuffer, numberOfFiles());
    
    //Create Nx * Ny * Nz array of zero's, this array will be the kernel array denoted Dk.
    float * DaI = malloc(stackSize);
    
    for (int i = 0; i < getImageSize()*numberOfFiles(); i++) {
        DaI[i] = 0;
    }
    
    //Apply the 3D FFT to Da.
        //i) Apply the 1-D FFT to each row of Da, denote this result Da1.
    
    // ------------------------> FFT row-wise <------------------------
    float * Da1R = malloc(stackSize);   //Da1R real component,
    float * Da1I = malloc(stackSize);   //Da1I imaginary component

    for (int i = 0; i < getImageHeight()*numberOfFiles(); i++) {
        int offset = i * getImageWidth();


        //memory copy row of data into Da1
        memcpy(Da1R+offset, DaR+offset, getImageWidth());
        memcpy(Da1I+offset, DaI+offset, getImageWidth());
        

        
        //transform in the forward direction row-wise
        FFT(1, ClosestPower2(getImageWidth()), Da1R+offset, Da1I+offset);
        //DFT(1, getImageWidth(), Da1R+offset, Da1I+offset);
        
        //printf("looped once!\n");

        //printf("Computed FFT at row number %i of image number %i\n", i%getImageWidth(), i%(getImageHeight()*numberOfFiles()));
    }
    printf("FFT row-wise done!\n");
    
    // ------------------------> FFT column-wise <------------------------ 
    float * Da2R = malloc(stackSize);   //Da2R real component,
    float * Da2I = malloc(stackSize);   //Da2I imaginary component
    
    for (int i = 0; i < getImageWidth()*numberOfFiles(); i++) {
        float* tmpR = malloc(getImageHeight()*sizeof(float));
        float* tmpI = malloc(getImageHeight()*sizeof(float));

        for (int j = 0; j < getImageHeight(); j++) {
            tmpR[j] = Da1R[i+(j*getImageWidth())];
            tmpI[j] = Da1I[i+(j*getImageWidth())];
        }
        
        FFT(FFT_FORWARD, ClosestPower2(getImageHeight()), tmpR, tmpI);
        //DFT(1, getImageHeight(), tmpR, tmpI);
        
        //set into new array
        for (int j = 0; j < getImageHeight(); j++) {
            Da2R[i+(j*getImageWidth())] = tmpR[j];
            Da2I[i+(j*getImageWidth())] = tmpI[j];
        }
    }
    printf("FFT column-wise done!\n");

    
    // ------------------------> FFT slice-wise <------------------------ 
    float * Da3R = malloc(stackSize);   //Da2R real component,
    float * Da3I = malloc(stackSize);   //Da2I imaginary component
    
    for (int i = 0; i < getImageWidth()*getImageHeight(); i++) {
        float* tmpR = malloc(numberOfFiles()*sizeof(float));
        float* tmpI = malloc(numberOfFiles()*sizeof(float));
        
        for (int j = 0; j < numberOfFiles(); j++) {
            tmpR[j] = Da2R[i+(j*getImageWidth()*getImageHeight())];
            tmpI[j] = Da2I[i+(j*getImageWidth()*getImageHeight())];
        }
        //NewFFT(8, tmpR, tmpI);

        DFT(1, numberOfFiles(), tmpR, tmpI);
        
        //set into new array
        for (int j = 0; j < numberOfFiles(); j++) {
            Da3R[i+(j*getImageWidth()*getImageHeight())] = tmpR[j];
            Da3I[i+(j*getImageWidth()*getImageHeight())] = tmpI[j];
        }
    }
    printf("FFT slice-wise done!\n");

    // ------------------------> Divide Da3 by (Nx*Ny*Nz) denoted DaFFT3d <------------------------ 
    float * DaFFT3dR = malloc(stackSize);
    float * DaFFT3dI = malloc(stackSize);
    
    for (int i = 0; i < getImageHeight()*getImageWidth()*numberOfFiles(); i++) {
        DaFFT3dR[i] = Da3R[i]/(getImageHeight()*getImageWidth()*numberOfFiles());
        DaFFT3dI[i] = Da3I[i]/(getImageHeight()*getImageWidth()*numberOfFiles());
    }
    
    // ------------------------> Build kernel denoted Dk <----------------------------------
    float* DkR = malloc(stackSize);
    float* DkI = malloc(stackSize);

    //generate Laplacian
    for (int i = 0; i < numberOfFiles()*getImageWidth()*getImageHeight(); i++) {
        if (i == ((numberOfFiles()*getImageWidth()*getImageHeight())/2)) {
            DkR[i] = 0;
        }
        else{
            DkR[i] = -1;
        }
        DkI[i] = 0;
    }
    
    //    float filtX[3] = {-1, 0, 1};
    //    float filtY[3] = {-1, 0, 1};
    //    float filtZ[3] = {-1, 0, 1};
    //    
    //    float dvfXYZ[3][3][3];
    //    float dvfYZX[3][3][3];
    //    float dvfZXY[3][3][3];
    //    
    //    for(int i = 0; i < 3; i++){
    //        for(int j = 0; j < 3; j++){
    //            for(int k = 0; k < 3; k++){
    //                dvfXYZ[i][j][k] = - filtX[i] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
    //                dvfYZX[i][j][k] = - filtY[j] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
    //                dvfZXY[i][j][k] = - filtZ[k] * exp(-((pow(filtX[i],2)+pow(filtY[j],2)+pow(filtZ[k],2))/2));
    //            }
    //        }
    //    }
    
    
// FFT the kernel now
    
    // ------------------------> FFT row-wise <------------------------
    float * Dk1R = malloc(stackSize);   //Da1R real component,
    float * Dk1I = malloc(stackSize);   //Da1I imaginary component
    
    for (int i = 0; i < getImageHeight()*numberOfFiles(); i++) {
        int offset = i * getImageWidth();
        
        //memory copy row of data into Da1
        memcpy(Dk1R+offset, DkR+offset, getImageWidth());
        memcpy(Dk1I+offset, DkI+offset, getImageWidth());
        
        //transform in the forward direction row-wise
        //FFT(FFT_FORWARD, ClosestPower2(getImageWidth()), Da1R+offset, Da1I+offset);
        DFT(1, getImageWidth(), Dk1R+offset, Dk1I+offset);

        //printf("looped once!\n");
        
        //printf("Computed FFT at row number %i of image number %i\n", i%getImageWidth(), i%(getImageHeight()*numberOfFiles()));
    }
    printf("kernel FFT row-wise done!\n");
    
    
    
    // ------------------------> FFT column-wise <------------------------ 
    float * Dk2R = malloc(stackSize);   //Da2R real component,
    float * Dk2I = malloc(stackSize);   //Da2I imaginary component
    
    for (int i = 0; i < getImageWidth()*numberOfFiles(); i++) {
        float* tmpR = malloc(getImageHeight()*sizeof(float));
        float* tmpI = malloc(getImageHeight()*sizeof(float));
        
        for (int j = 0; j < getImageHeight(); j++) {
            tmpR[j] = Dk1R[i+(j*getImageWidth())];
            tmpI[j] = Dk1I[i+(j*getImageWidth())];
        }

        //FFT(FFT_FORWARD, ClosestPower2(getImageHeight()), tmpR, tmpI);
        DFT(1, getImageHeight(), tmpR, tmpI);

        //set into new array
        for (int j = 0; j < getImageHeight(); j++) {
            Dk2R[i+(j*getImageWidth())] = tmpR[j];
            Dk2I[i+(j*getImageWidth())] = tmpI[j];
        }
    }
    printf("kernel FFT column-wise done!\n");
    
    
    // ------------------------> FFT slice-wise <------------------------ 
    float * Dk3R = malloc(stackSize);   //Da2R real component,
    float * Dk3I = malloc(stackSize);   //Da2I imaginary component
    
    for (int i = 0; i < getImageWidth()*getImageHeight(); i++) {
        float* tmpR = malloc(numberOfFiles()*sizeof(float));
        float* tmpI = malloc(numberOfFiles()*sizeof(float));
        
        for (int j = 0; j < numberOfFiles(); j++) {
            tmpR[j] = Dk2R[i+(j*getImageWidth()*getImageHeight())];
            tmpI[j] = Dk2I[i+(j*getImageWidth()*getImageHeight())];
        }
        
        DFT(1, numberOfFiles(), tmpR, tmpI);
        
        //set into new array
        for (int j = 0; j < numberOfFiles(); j++) {
            Dk3R[i+(j*getImageWidth()*getImageHeight())] = tmpR[j];
            Dk3I[i+(j*getImageWidth()*getImageHeight())] = tmpI[j];
        }
    }
    printf("kernel FFT slice-wise done!\n");
    
//kernel FFT done!
    
    // ------------------------> Divide Dk3 by (Nx*Ny*Nz) denoted DkFFT3d <------------------------ 
    float * DkFFT3dR = malloc(stackSize);
    float * DkFFT3dI = malloc(stackSize);
    
    for (int i = 0; i < getImageHeight()*getImageWidth()*numberOfFiles(); i++) {
        DkFFT3dR[i] = Dk3R[i]/(getImageHeight()*getImageWidth()*numberOfFiles());
        DkFFT3dI[i] = Dk3I[i]/(getImageHeight()*getImageWidth()*numberOfFiles());
    }
    
    // ------------------------> Take the complex conjugate of DaFFT3d <---------------------------
    for(int i = 0; i < getImageWidth()*getImageHeight()*numberOfFiles(); i++){
        DaFFT3dI[i] = -DaFFT3dI[i];
    }
    
    // ------------------------> (Convolution) Multiply DaFFT3d conjugate by DkFFT3d <-------------
    for(int i = 0; i < getImageWidth()*getImageHeight()*numberOfFiles(); i++){
//        DaFFT3dR[i] = DaFFT3dR[i] * DkFFT3dR[i];
//        DaFFT3dI[i] = DaFFT3dI[i] * DkFFT3dI[i];
        
        //due to a trick of the tail, store the results in Da3 instead
        Da3R[i] = DaFFT3dR[i] * DkFFT3dR[i];
        Da3I[i] = DaFFT3dI[i] * DkFFT3dI[i];
    }
    printf("kernel Convolution done!\n");
    
    // ------------------------> Inverse FFT slice-wise <------------------------ 
    for (int i = 0; i < getImageWidth()*getImageHeight(); i++) {
        float* tmpR = malloc(numberOfFiles()*sizeof(float));
        float* tmpI = malloc(numberOfFiles()*sizeof(float));
        
        for (int j = 0; j < numberOfFiles(); j++) {
            tmpR[j] = Da3R[i+(j*getImageWidth()*getImageHeight())];
            tmpI[j] = Da3I[i+(j*getImageWidth()*getImageHeight())];
        }
        //NewFFT(8, tmpR, tmpI);

        DFT(-1, numberOfFiles(), tmpR, tmpI);
        
        //set into new array
        for (int j = 0; j < numberOfFiles(); j++) {
            Da2R[i+(j*getImageWidth()*getImageHeight())] = tmpR[j];
            Da2I[i+(j*getImageWidth()*getImageHeight())] = tmpI[j];
        }
    }
    printf("IFFT slice-wise done!\n");

    // ------------------------> Inverse FFT column-wise <------------------------ 
    
    for (int i = 0; i < getImageWidth()*numberOfFiles(); i++) {
        float* tmpR = malloc(getImageHeight()*sizeof(float));
        float* tmpI = malloc(getImageHeight()*sizeof(float));
        
        for (int j = 0; j < getImageHeight(); j++) {
            tmpR[j] = Da2R[i+(j*getImageWidth())];
            tmpI[j] = Da2I[i+(j*getImageWidth())];
        }
        
        FFT(FFT_REVERSE, ClosestPower2(getImageHeight()), tmpR, tmpI);

        //DFT(-1, getImageHeight(), tmpR, tmpI);
        
        //set into new array
        for (int j = 0; j < getImageHeight(); j++) {
            Da1R[i+(j*getImageWidth())] = tmpR[j];
            Da1I[i+(j*getImageWidth())] = tmpI[j];
        }
    }
    printf("IFFT column-wise done!\n");
    
    // ------------------------> Inverse FFT row-wise <------------------------ 

    for (int i = 0; i < getImageHeight()*numberOfFiles(); i++) {
        int offset = i * getImageWidth();
        
        //transform in the forward direction row-wise
        FFT(-1, 8, Da1R+offset, Da1I+offset);
        //DFT(-1, getImageWidth(), Da1R+offset, Da1I+offset);
        //NewFFT(8, Da1R+offset, Da1I+offset);

        //memory copy row of data into Da1
        memcpy(DaR+offset, Da1R+offset, getImageWidth());
        memcpy(DaI+offset, Da1I+offset, getImageWidth());
        
        //printf("looped once!\n");
        
        //printf("Computed FFT at row number %i of image number %i\n", i%getImageWidth(), i%(getImageHeight()*numberOfFiles()));
    }
    printf("IFFT row-wise done!\n");

    

    //collect result ready for writing
    bigBuffer = denormaliseStack(DaR, numberOfFiles());
    

    /* ----------------------------------------> KERNEL ENDS HERE! <---------------------------------------- */
    
    //save all images from buffer
    for (int i = 0; i < depth; i++) {
        uint8 *buffer = malloc(sizeof(uint8)*getImageSize());
        memcpy(buffer, bigBuffer+(i*getImageSize()), getImageSize());
                
        char* file = substring(locateLast(outputImageFileName,'/')+1, strlen(outputImageFileName), outputImageFileName);

        //printf("filename : %s", file);
        
        char* path = substring(0, locateLast(outputImageFileName,'/')+1, outputImageFileName);
        
        char* cutDownFile = substring(0, locateLast(file,'.'), file);
        char* extension = substring(locateLast(file,'.'), (int)strlen(file),file);

        char* newName = cutDownFile;
        char numericalRepresentation[200];
        sprintf(numericalRepresentation, "%d", i);
        newName = strcat(newName, numericalRepresentation);
        newName = strcat(newName, extension);
        newName = strcat(path, newName);
        
        saveImage(newName, buffer);   
    }
    
#else
    
//#define DebugPrintArrayValues // <- uncomment this if you want to print arrays to verify transform results
        //timer variables
    clock_t startTime, stopTime;

    int textpow = 6;
    
        //initialize and set test arrays
    float* slowDataReal = malloc(sizeof(float)*pow(2, textpow));
    float* slowDataImag = malloc(sizeof(float)*pow(2, textpow));
    
    for (int i = 0; i < pow(2, textpow); i++) {
        slowDataReal[i] = i;
        slowDataImag[i] = i;
    }
    
    float* dataReal = malloc(sizeof(float)*pow(2, textpow));
    float* dataImag = malloc(sizeof(float)*pow(2, textpow));
    
    for (int i = 0; i < pow(2, textpow); i++) {
        dataReal[i] = i;
        dataImag[i] = i;
    }
    
    
#ifdef DebugPrintArrayValues
        //verify that all four arrays populated correctly
    for (int i = 0; i < pow(2, textpow); i++) {
        printf("before dft: real data at index %i contains value \t\t:%f\n", i, slowDataReal[i]);
        printf("before dft: imaginary data at index %i contains value \t:%f\n", i, slowDataImag[i]);
        printf("\n");
        printf("before fft: real data at index %i contains value \t\t:%f\n", i, slowDataReal[i]);
        printf("before fft: imaginary data at index %i contains value \t:%f\n", i, slowDataImag[i]);
        printf("\n");
        printf("\n");
    }
#endif
    
    //time forward DFT execution
    startTime = clock();
    DFT(1, pow(2, textpow), slowDataReal, slowDataImag);
    stopTime = clock();
    printf("Time to perform DFT was %f seconds\n", (double)(stopTime-startTime)/CLOCKS_PER_SEC);
    
    //time forward FFT execution
    startTime = clock();
    FFT(FFT_FORWARD,  textpow, dataReal, dataImag);
    stopTime = clock();
    printf("Time to perform FFT was %f seconds\n", (double)(stopTime-startTime)/CLOCKS_PER_SEC);

#ifdef DebugPrintArrayValues
    //verify that all four arrays have been transformed
    for (int i = 0; i < pow(2, textpow); i++) {
        printf("after dft: real data at index %i contains value \t\t:%f\n", i, slowDataReal[i]);
        printf("after dft: imaginary data at index %i contains value \t:%f\n", i, slowDataImag[i]);
        printf("\n");
        printf("after fft: real data at index %i contains value \t\t:%f\n", i, dataReal[i]);
        printf("after fft: imaginary data at index %i contains value \t:%f\n", i, dataImag[i]);
        printf("\n");
        printf("\n");
    }
#endif
    
    //time reverse DFT execution
    startTime = clock();
    DFT(-1, pow(2, textpow), slowDataReal, slowDataImag);
    stopTime = clock();
    printf("Time to perform Inverse DFT was %f seconds\n", (double)(stopTime-startTime)/CLOCKS_PER_SEC);
    
    //time reverse FFT execution
    startTime = clock();
    FFT(FFT_REVERSE, textpow, dataReal, dataImag);
    stopTime = clock();
    printf("Time to perform Inverse FFT was %f seconds\n", (double)(stopTime-startTime)/CLOCKS_PER_SEC);
    
#ifdef DebugPrintArrayValues
    //verify that all four arrays have inverse transform corresponding closely to their original values.
    for (int i = 0; i < pow(2, textpow); i++) {
        printf("after inv dft: real data at index %i contains value \t\t:%f\n", i, slowDataReal[i]);
        printf("after inv dft: imaginary data at index %i contains value \t:%f\n", i, slowDataImag[i]);
        printf("\n");
        printf("after inv fft: real data at index %i contains value \t\t:%f\n", i, dataReal[i]);
        printf("after inv fft: imaginary data at index %i contains value \t:%f\n", i, dataImag[i]);
        printf("\n");
        printf("\n");
    }
#endif
    
    
#endif
    return 0;
}

