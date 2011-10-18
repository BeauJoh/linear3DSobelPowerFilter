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
    
    
#ifdef complexArrayTesting
    /* three dimensional array testing */
    
    float testit[4][4][4];
    float testitI[4][4][4];
    
    float z = 0;
    for (int i = 0; i < 4; i ++) {
        for (int j = 0; j < 4; j ++) {
            for (int k = 0; k < 4; k ++) {
                
                if (k == 3 || j == 3 || i == 3) {
                    testit[i][j][k] = 0;
                    testitI[i][j][k] = 0;
                }
                else{
                    testit[i][j][k] = z;
                    testitI[i][j][k] = 0;
                    z++;
                }
                
                
            }
        }
    }
    
    printf("Before convolution");
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                printf("%f \n", testit[i][j][k]);
            }
        }
    }

    //row wise fft
    for (int i = 0; i < 4; i ++) {
        for (int j = 0; j < 4; j ++) {
            float * tmpRowR = malloc(sizeof(float)*4);
            float * tmpRowI = malloc(sizeof(float)*4);
            
            //collect a row
            for (int k = 0; k < 4; k ++) {
                // throw into a tmp array to do FFT upon
                tmpRowR[k] = testit[i][j][k];
                tmpRowI[k] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_FORWARD, 2, tmpRowR, tmpRowI);
            
            // store the resulting row into original array
            for (int k = 0; k < 4; k ++) {
                // throw into a tmp array to do FFT upon
                testit[i][j][k] = tmpRowR[k];
                testitI[i][j][k] = tmpRowI[k];
            }
        }
    }
    
    //column wise fft
    for (int i = 0; i < 4; i ++) {
        for (int k = 0; k < 4; k ++) {
            float * tmpColR = malloc(sizeof(float)*4);
            float * tmpColI = malloc(sizeof(float)*4);
            
            for (int j = 0; j < 4; j ++) {
                // throw into a tmp array to do FFT upon
                tmpColR[j] = testit[i][j][k];
                tmpColI[j] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_FORWARD, 2, tmpColR, tmpColI);
            
            for (int j = 0; j < 4; j ++) {
                // throw into a tmp array to do FFT upon
                testit[i][j][k] = tmpColR[j];
                testitI[i][j][k] = tmpColI[j];
            }
        }
    }
    
    //slice wise fft
    for (int j = 0; j < 4; j ++) {
        for (int k = 0; k < 4; k ++) {
            float * tmpSliR = malloc(sizeof(float)*4);
            float * tmpSliI = malloc(sizeof(float)*4);
            
            //throw present slice into a tmp array to do FFT upon
            for (int i = 0; i < 4; i ++) {
                tmpSliR[i] = testit[i][j][k];
                tmpSliI[i] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_FORWARD, 2, tmpSliR, tmpSliI);
            
            //collect present slice into original 4*4*4 array
            for (int i = 0; i < 4; i ++) {
                testit[i][j][k] = tmpSliR[i];
                testitI[i][j][k] = tmpSliI[i];
            }
            
        }
    }
    
    
    //slice wise inverse fft
    for (int j = 0; j < 4; j ++) {
        for (int k = 0; k < 4; k ++) {
            float * tmpSliR = malloc(sizeof(float)*4);
            float * tmpSliI = malloc(sizeof(float)*4);
            
            //throw present slice into a tmp array to do FFT upon
            for (int i = 0; i < 4; i ++) {
                tmpSliR[i] = testit[i][j][k];
                tmpSliI[i] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_REVERSE, 2, tmpSliR, tmpSliI);
            
            //collect present slice into original 4*4*4 array
            for (int i = 0; i < 4; i ++) {
                testit[i][j][k] = tmpSliR[i];
                testitI[i][j][k] = tmpSliI[i];
            }
            
        }
    }
    
    //column wise inverse fft
    for (int i = 0; i < 4; i ++) {
        for (int k = 0; k < 4; k ++) {
            float * tmpColR = malloc(sizeof(float)*4);
            float * tmpColI = malloc(sizeof(float)*4);
            
            for (int j = 0; j < 4; j ++) {
                // throw into a tmp array to do FFT upon
                tmpColR[j] = testit[i][j][k];
                tmpColI[j] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_REVERSE, 2, tmpColR, tmpColI);
            
            for (int j = 0; j < 4; j ++) {
                // throw into a tmp array to do FFT upon
                testit[i][j][k] = tmpColR[j];
                testitI[i][j][k] = tmpColI[j];
            }
        }
    }

    
    printf("After convolution");
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                printf("%f \n", testit[i][j][k]);
            }
        }
    }

    //row wise inv fft
    for (int i = 0; i < 4; i ++) {
        for (int j = 0; j < 4; j ++) {
            float * tmpRowR = malloc(sizeof(float)*4);
            float * tmpRowI = malloc(sizeof(float)*4);
            
            //collect a row
            for (int k = 0; k < 4; k ++) {
                // throw into a tmp array to do FFT upon
                tmpRowR[k] = testit[i][j][k];
                tmpRowI[k] = testitI[i][j][k];
            }
            
            //apply FFT
            FFT(FFT_REVERSE, 2, tmpRowR, tmpRowI);
            
            // store the resulting row into original array
            for (int k = 0; k < 4; k ++) {
                // throw into a tmp array to do FFT upon
                testit[i][j][k] = tmpRowR[k];
                testitI[i][j][k] = tmpRowI[k];
            }
        }
    }
    
    
    printf("After inverse convolution");
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                printf("%f \n", testit[i][j][k]);
            }
        }
    }
    
    /*
    printf("row-wise traversal \n");
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                printf("at index [%i][%i][%i] -> %f \n", i, j, k, testit[i][j][k]);
            }
        }
    }
    
    printf("column-wise traversal \n");
    for (int i = 0; i < 3; i ++) {
        for (int k = 0; k < 3; k ++) {
            for (int j = 0; j < 3; j ++) {
                printf("at index [%i][%i][%i] -> %f \n", i, j, k, testit[i][j][k]);
            }
        }
    }
    
    printf("slice-wise traversal \n");
    for (int j = 0; j < 3; j ++) {
        for (int k = 0; k < 3; k ++) {
            for (int i = 0; i < 3; i ++) {
                printf("at index [%i][%i][%i] -> %f \n", i, j, k, testit[i][j][k]);
            }
        }
    }
    */
#endif
    
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
    
    //time it 
    clock_t startTime, stopTime;
    startTime = clock();

    //pull down 3*3*3 tiles from the entire image
    for (int z = 0; z < numberOfFiles()-3; z +=3) {
        for (int y = 0; y < getImageHeight(); y +=3) {
            for (int x = 0; x < getImageWidth(); x +=3) {
              for (int c = 0; c < getSamplesPerPixel()-1; c++) {
                
                    float DatR[4][4][4];
                    float DatI[4][4][4];
                    
                    
                    //set out 4 window size (need 3 for correct filter focus, however need to be dyadic for fft, solution use window 
                    //size of 4*4*4 whilst only populating the 3*3*3 and padding the rest with zeros)
                    for(int i = 0; i < 4; i++){
                        for(int j = 0; j < 4; j++){
                            for(int k = 0; k < 4; k++){
                                if (k == 3 || j == 3 || i == 3) {
                                    DatR[i][j][k] = 0;
                                    DatI[i][j][k] = 0;
                                }
                                else{
                                    DatR[i][j][k] = DaR[((z+i)*getImageHeight()*getImageWidth()*getSamplesPerPixel()) + ((y+j)*getImageWidth()*getSamplesPerPixel()) + (x+k)*getSamplesPerPixel()+c];
                                    DatI[i][j][k] = DaI[((z+i)*getImageHeight()*getImageWidth()*getSamplesPerPixel()) + ((y+j)*getImageWidth()*getSamplesPerPixel()) + (x+k)*getSamplesPerPixel()+c];
                                }
                            }
                        }
                    }
                    
                    
                    //row wise fft
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            float * tmpRowR = malloc(sizeof(float)*4);
                            float * tmpRowI = malloc(sizeof(float)*4);
                            
                            //collect a row
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                tmpRowR[k] = DatR[i][j][k];
                                tmpRowI[k] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpRowR, tmpRowI);
                            
                            // store the resulting row into original array
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                DatR[i][j][k] = tmpRowR[k];
                                DatI[i][j][k] = tmpRowI[k];
                            }
                            free(tmpRowR);
                            free(tmpRowI);
                        }
                    }
                    
                    //column wise fft
                    for (int i = 0; i < 4; i ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpColR = malloc(sizeof(float)*4);
                            float * tmpColI = malloc(sizeof(float)*4);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                tmpColR[j] = DatR[i][j][k];
                                tmpColI[j] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpColR, tmpColI);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                DatR[i][j][k] = tmpColR[j];
                                DatI[i][j][k] = tmpColI[j];
                            }
                            free(tmpColR);
                            free(tmpColI);
                        }
                    }
                    
                    //slice wise fft
                    for (int j = 0; j < 4; j ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpSliR = malloc(sizeof(float)*4);
                            float * tmpSliI = malloc(sizeof(float)*4);
                            
                            //throw present slice into a tmp array to do FFT upon
                            for (int i = 0; i < 4; i ++) {
                                tmpSliR[i] = DatR[i][j][k];
                                tmpSliI[i] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpSliR, tmpSliI);
                            
                            //collect present slice into original 4*4*4 array
                            for (int i = 0; i < 4; i ++) {
                                DatR[i][j][k] = tmpSliR[i];
                                DatI[i][j][k] = tmpSliI[i];
                            }
                            
                            free(tmpSliR);
                            free(tmpSliI);
                        }
                    }

                    
                    // ------------------------> Divide Da by (3*3*3) denoted Dk <------------------------ 
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            for (int k = 0; k < 4; k ++) {
                                DatR[i][j][k] = DatR[i][j][k] / (4*4*4);
                                DatI[i][j][k] = DatI[i][j][k] / (4*4*4);
                            }
                        }
                    }

            
                    
                    //convolution 
                    //generate the kernel
                    //(Sobel Power Filter Bank)
                    float DkR[4][4][4];
                    float DkI[4][4][4];
                    
                    float filtX[3] = {-1, 0, 1};
                    float filtY[3] = {-1, 0, 1};
                    float filtZ[3] = {-1, 0, 1};
                    
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            for (int k = 0; k < 4; k++) {
                                if (i == 3 || j == 3 || k == 3) {
                                    DkR[i][j][k] = 0;
                                }
                                else {
                                    if ((k == 0) && (i == 0) && (j == 0)){
                                        DkR[i][j][k] = 0.223;
                                    }
                                    else if ((i == 0) && (j == 2)){
                                        DkR[i][j][k] = 0.223;
                                    }
                                    else if ((i == 2) && (j == 0)){
                                        DkR[i][j][k] = -0.223;
                                    }
                                    else if ((i == 2) && (j == 2)){
                                        DkR[i][j][k] = -0.223;
                                    }
                                    else{
                                        DkR[i][j][k] = 0;
                                    }
                                     
                                }
                                
                                //printf("%f\n", DkR[i][j][k]);
                                DkI[i][j][k] = 0;
                                
                            }
                        }
                    }
                    
                    //Apply forward transform upon filter
                    //First x-wise
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            float * tmpRowR = malloc(sizeof(float)*4);
                            float * tmpRowI = malloc(sizeof(float)*4);
                            
                            //collect a row
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                tmpRowR[k] = DkR[i][j][k];
                                tmpRowI[k] = DkI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpRowR, tmpRowI);
                            
                            // store the resulting row into original array
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                DkR[i][j][k] = tmpRowR[k];
                                DkI[i][j][k] = tmpRowI[k];
                            }
                            
                            free(tmpRowR);
                            free(tmpRowI);
                        }
                    }
                    
                    //Then y-wise
                    for (int i = 0; i < 4; i ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpColR = malloc(sizeof(float)*4);
                            float * tmpColI = malloc(sizeof(float)*4);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                tmpColR[j] = DkR[i][j][k];
                                tmpColI[j] = DkI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpColR, tmpColI);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                DkR[i][j][k] = tmpColR[j];
                                DkI[i][j][k] = tmpColI[j];
                            }
                            
                            free(tmpColR);
                            free(tmpColI);
                        }
                    }
                    
                    //Then z-wise
                    for (int j = 0; j < 4; j ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpSliR = malloc(sizeof(float)*4);
                            float * tmpSliI = malloc(sizeof(float)*4);
                            
                            //throw present slice into a tmp array to do FFT upon
                            for (int i = 0; i < 4; i ++) {
                                tmpSliR[i] = DkR[i][j][k];
                                tmpSliI[i] = DkI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_FORWARD, 2, tmpSliR, tmpSliI);
                            
                            //collect present slice into original 4*4*4 array
                            for (int i = 0; i < 4; i ++) {
                                DkR[i][j][k] = tmpSliR[i];
                                DkI[i][j][k] = tmpSliI[i];
                            }
                            free(tmpSliR);
                            free(tmpSliI);
                        }
                    }
                    
                    if (z==0&&y==0&&x==0) {
                        for (int i = 0; i < 3; i ++) {
                            for (int j = 0; j < 3; j ++) {
                                for (int k = 0; k < 3; k ++) {
                                    printf("at index [%i][%i][%i] -> %f \n", i, j, k, DkR[i][j][k]);
                                }
                            }
                        }
                        
                    }

                    
                    //apply convolution
                    // ------------------------> Divide Dk by (3*3*3) denoted Dk <------------------------ 
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            for (int k = 0; k < 4; k ++) {
                                DkR[i][j][k] = DkR[i][j][k] / (4*4*4);
                                DkI[i][j][k] = DkI[i][j][k] / (4*4*4);
                            }
                        }
                    }
                    
                    // ------------------------> Take the complex conjugate of Da <---------------------------
                    for (int i = 0; i < 3; i ++) {
                        for (int j = 0; j < 3; j ++) {
                            for (int k = 0; k < 3; k ++) {
                                DatI[i][j][k] = -DatI[i][j][k];
                            }
                        }
                    }

                    // ------------------------> (Convolution) Multiply Da conjugate by Dk <-------------
                    for (int i = 0; i < 3; i ++) {
                        for (int j = 0; j < 3; j ++) {
                            for (int k = 0; k < 3; k ++) {
                                DatR[i][j][k] = DatR[i][j][k] * DkR[i][j][k];
                                DatI[i][j][k] = DatI[i][j][k] * DkI[i][j][k];
                            }
                        }
                    }
                    //end of convolution
                    
                    //inverse transformation
                    //First z-wise (slice)
                    for (int j = 0; j < 4; j ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpSliR = malloc(sizeof(float)*4);
                            float * tmpSliI = malloc(sizeof(float)*4);
                            
                            //throw present slice into a tmp array to do FFT upon
                            for (int i = 0; i < 4; i ++) {
                                tmpSliR[i] = DatR[i][j][k];
                                tmpSliI[i] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_REVERSE, 2, tmpSliR, tmpSliI);
                            
                            //collect present slice into original 4*4*4 array
                            for (int i = 0; i < 4; i ++) {
                                DatR[i][j][k] = tmpSliR[i];
                                DatI[i][j][k] = tmpSliI[i];
                            }
                            free(tmpSliR);
                            free(tmpSliI);
                        }
                    }
                    
                    //Then y-wise (column wise inverse fft)
                    for (int i = 0; i < 4; i ++) {
                        for (int k = 0; k < 4; k ++) {
                            float * tmpColR = malloc(sizeof(float)*4);
                            float * tmpColI = malloc(sizeof(float)*4);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                tmpColR[j] = DatR[i][j][k];
                                tmpColI[j] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_REVERSE, 2, tmpColR, tmpColI);
                            
                            for (int j = 0; j < 4; j ++) {
                                // throw into a tmp array to do FFT upon
                                DatR[i][j][k] = tmpColR[j];
                                DatI[i][j][k] = tmpColI[j];
                            }
                            free(tmpColR);
                            free(tmpColI);
                        }
                    }

                    
                    //Finally x-wise (row wise inv fft)
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            float * tmpRowR = malloc(sizeof(float)*4);
                            float * tmpRowI = malloc(sizeof(float)*4);
                            
                            //collect a row
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                tmpRowR[k] = DatR[i][j][k];
                                tmpRowI[k] = DatI[i][j][k];
                            }
                            
                            //apply FFT
                            FFT(FFT_REVERSE, 2, tmpRowR, tmpRowI);
                            
                            // store the resulting row into original array
                            for (int k = 0; k < 4; k ++) {
                                // throw into a tmp array to do FFT upon
                                DatR[i][j][k] = tmpRowR[k];
                                DatI[i][j][k] = tmpRowI[k];
                         }
                            free(tmpRowR);
                            free(tmpRowI);
                        }
                    }                
                    
                    // ------------------------> Multiply Da by (3*3*3) denoted Da' <------------------------ 
                    for (int i = 0; i < 4; i ++) {
                        for (int j = 0; j < 4; j ++) {
                            for (int k = 0; k < 4; k ++) {
                                DatR[i][j][k] = DatR[i][j][k] * (4*4*4);
                                DatI[i][j][k] = DatI[i][j][k] * (4*4*4);
                            }
                        }
                    }
                    
                    //finally store our results back into original DaR array
                    for(int i = 0; i < 3; i++){
                        for(int j = 0; j < 3; j++){
                            for(int k = 0; k < 3; k++){
                                DaR[((z+i)*getImageHeight()*getImageWidth()*getSamplesPerPixel()) + ((y+j)*getImageWidth()*getSamplesPerPixel()) + (x+k)*getSamplesPerPixel()+c] = DatR[i][j][k];
                                DaI[((z+i)*getImageHeight()*getImageWidth()*getSamplesPerPixel()) + ((y+j)*getImageWidth()*getSamplesPerPixel()) + (x+k)*getSamplesPerPixel()+c] = DatI[i][j][k];                            
                        }
                    }
                }
            }
        }
    }
    }
    
    // stop timer and show times
    stopTime = clock();
    printf("Time to perform convolution was %f seconds\n", (double)(stopTime-startTime)/CLOCKS_PER_SEC);
    
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
    

    return 0;
}

