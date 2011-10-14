/*
 *  RGBAUtilities.c
 *  MedialAxisTransform
 *
 *
 *  Created by Beau Johnston on 13/07/11.
 *  Copyright (C) 2011 by Beau Johnston.
 *
 *  Please email me if you have any comments, suggestions or advice:
 *                              beau@inbeta.org
 *
 *  read_png_file & write_png_file functionality, Guillaume Cottenceau, 2002-2010
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


#include "RGBAUtilities.h"

int x, y;
png_structp pngPtr;
png_infop infoPtr;
int numberOfPasses;
png_bytep * rowPointers;
uint32 imageLength, imageWidth, config, bitsPerSample, samplesPerPixel, bitsPerPixel, imageBitSize, imageSize;
uint64 linebytes;

union FloatAndByte
{
    float   f;
    uint8   c[0];
};


void abort_(const char * s, ...)
{
    va_list args;
    va_start(args, s);
    vfprintf(stderr, s, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();
}


void readPngFile(char* file_name)
{
    //printf("%s\n", file_name);
    char header[8];    // 8 is the maximum size that can be checked
    
    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp){
        printf("Unable to open file '%s'\n", file_name);
        exit(-1);
        //abort_("[read_png_file] File %s could not be opened for reading", file_name);
    }
    fread(header, 1, 8, fp);
//    if (png_sig_cmp(header, 0, 8))
//        abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);
//    
    
    /* initialize stuff */
    pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!pngPtr)
        abort_("[read_png_file] png_create_read_struct failed");
    
    infoPtr = png_create_info_struct(pngPtr);
    if (!infoPtr)
        abort_("[read_png_file] png_create_info_struct failed");
    
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[read_png_file] Error during init_io");
    
    png_init_io(pngPtr, fp);
    png_set_sig_bytes(pngPtr, 8);
    
    png_read_info(pngPtr, infoPtr);
    
    imageWidth = png_get_image_width(pngPtr, infoPtr);
    imageLength = png_get_image_height(pngPtr, infoPtr);
    config = png_get_color_type(pngPtr, infoPtr);
    bitsPerSample = png_get_bit_depth(pngPtr, infoPtr); // = 8 bits
        
    numberOfPasses = png_set_interlace_handling(pngPtr);
    png_read_update_info(pngPtr, infoPtr);
    
    imageWidth = png_get_image_width(pngPtr, infoPtr);
    imageLength = png_get_image_height(pngPtr, infoPtr);

    samplesPerPixel = png_get_channels(pngPtr, infoPtr); // = 4 bytes

    bitsPerPixel = samplesPerPixel*bitsPerSample;
    linebytes = samplesPerPixel * imageWidth; // = 640
    //linebytes = png_get_rowbytes(pngPtr, infoPtr); = 640
    imageBitSize = (sizeof(uint8) * imageWidth * imageLength * samplesPerPixel);
    imageSize = imageWidth * imageLength * samplesPerPixel;
    //printf("linebytes = %i, expected %i\n",linebytes,png_get_rowbytes(pngPtr, infoPtr));
    //printf("Image Height is %d", sizeof(png_bytep) * imageLength);

    
    /* read file */
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[read_png_file] Error during read_image");
    
    rowPointers = (png_bytep*) malloc(sizeof(png_bytep) * imageLength);
    for (y=0; y<imageLength; y++)
        rowPointers[y] = (png_byte*) malloc(png_get_rowbytes(pngPtr,infoPtr));
    
    png_read_image(pngPtr, rowPointers);
    
    fclose(fp);
}


void writePngFile(char* file_name)
{
    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        abort_("[write_png_file] File %s could not be opened for writing", file_name);
    
    
    /* initialize stuff */
    pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!pngPtr)
        abort_("[write_png_file] png_create_write_struct failed");
    
    infoPtr = png_create_info_struct(pngPtr);
    if (!infoPtr)
        abort_("[write_png_file] png_create_info_struct failed");
    
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[write_png_file] Error during init_io");
    
    png_init_io(pngPtr, fp);
    
    
    /* write header */
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[write_png_file] Error during writing header");
    
    png_set_IHDR(pngPtr, infoPtr, imageWidth, imageLength,
                 bitsPerSample, config, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
    png_write_info(pngPtr, infoPtr);
    
    
    /* write bytes */
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[write_png_file] Error during writing bytes");
    
    png_write_image(pngPtr, rowPointers);
    
    
    /* end write */
    if (setjmp(png_jmpbuf(pngPtr)))
        abort_("[write_png_file] Error during end of write");
    
    png_write_end(pngPtr, NULL);
        
    fclose(fp);
}

void cleanup(void){
    /* cleanup heap allocation */
    for (y=0; y<imageLength; y++)
        free(rowPointers[y]);
    free(rowPointers);
}

float* normalizeImage(uint8* input){
    //with 8 bits this obvously causes a rounding error, usually down to 0, solve this by storing as floats
    float* output = (float *)malloc(sizeof(float)*imageSize);
    
    for (int i = 0; i < imageSize; i++) {
        output[i] = ((float)input[i]/255.0f);
    }
    free(input);
    return output;
}

float* normaliseStack(uint8* input, int depth){
    //with 8 bits this obvously causes a rounding error, usually down to 0, solve this by storing as floats
    float* output = (float *)malloc(sizeof(float)*imageSize*depth);
    
    for (int i = 0; i < imageSize*depth; i++) {
        output[i] = ((float)input[i]/255.0f);
    }
    free(input);
    return output;
}

uint8* denormaliseStack(float* input, int depth){
    //with 8 bits this obvously causes a rounding error, usually down to 0, solve this by storing as floats
    uint8* output = (uint8 *)malloc(sizeof(uint8)*imageSize*depth);
    
    for (int i = 0; i < imageSize*depth; i++) {
        output[i] = ((float)input[i]*255.0f);
    }
    free(input);
    return output;
}

uint8* denormalizeImage(float*input){
    //with 8 bits this obvously causes a rounding error, usually down to 0, solve this by storing as floats
    uint8* output = (uint8 *)malloc(sizeof(uint8)*imageSize);
    for (int i = 0; i < imageSize; i++) {
        output[i] = ((uint8)input[i]*255.0f);
    }
    free(input);
    return output;
}

bool allPixelsAreNormal(uint8* image){
    #ifdef DEBUG
    printf("bits per sample is %i", bitsPerSample);
    #endif
    for (int i = 0; i < imageSize; i++) {
        if (image[i] > 1) {
            #ifdef DEBUG
            printf("Not normal at point %i with %d", i, image[i]);
            #endif
            return false;
        }
    }
    return true;
}

uint8* convolutedGetImage(void){
    uint8* image = (uint8*)malloc(sizeof(uint8) * imageLength);

    if (png_get_color_type(pngPtr, infoPtr) != PNG_COLOR_TYPE_RGBA)
        abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (%d) (is %d)",
                PNG_COLOR_TYPE_RGBA, png_get_color_type(pngPtr, infoPtr));
    
    for (y=0; y<imageLength; y++) {
        uint8* row = rowPointers[y];

        for (x=0; x<imageWidth; x++) {
            uint8* ptr = &(row[x*4]);
            image[y*imageWidth + x] = *ptr;
        }
    }
    
    return image;
}

uint8* getImage(void){
    uint8* image = (uint8*)malloc(sizeof(uint8) * imageWidth * imageLength * samplesPerPixel);
    for (y=0; y < imageLength; y++) {
        uint8* row = rowPointers[y];
        int origX = 0;
        for (x=0; x < linebytes; x+=4) {
            uint8* ptr = &(row[origX*4]);
            image[y*linebytes + x + 0] = ptr[0];
            image[y*linebytes + x + 1] = ptr[1];
            image[y*linebytes + x + 2] = ptr[2];
            image[y*linebytes + x + 3] = ptr[3];
            origX++;
        }
    }
    
    return image;
}

void convolutedSetImage(uint8* image){
    for (y=0; y<imageLength; y++) {
        uint8* row = rowPointers[y];
        for (x=0; x<imageWidth; x++) {
            uint8* ptr = &(row[x*4]);
            *ptr = image[y*imageWidth + x];
        }
    }
    
    free(image);
}

void clearImageBuffer(void){
    for (y=0; y<imageLength; y++) {
        uint8* row = rowPointers[y];
        int origX = 0;
        for (x=0; x<linebytes; x+=4) {
            uint8* ptr = &(row[origX*4]);
            ptr[0] = 0;
            ptr[1] = 0;
            ptr[2] = 0;
            ptr[3] = 0;
            origX++;
        }
    }
}

void setImage(uint8* image){
    
    for (y=0; y < imageLength; y++) {
        uint8* row = rowPointers[y];
        int origX = 0;
        for (x=0; x < linebytes; x+=4) {
            uint8* ptr = &(row[origX*4]);
            ptr[0] = image[y*linebytes + x + 0];
            ptr[1] = image[y*linebytes + x + 1];
            ptr[2] = image[y*linebytes + x + 2];
            ptr[3] = image[y*linebytes + x + 3];
            origX++;
        }
    }
    
    free(image);
}

float* norm(float* input, uint32 imageSize){
    float* output = malloc(sizeof(float)*imageSize);
    
    for (int i = 0; i < imageSize; i++) {
        output[i] = input[i]/255.0f;
    }
    return output;
}

float* denorm(float* input, uint32 imageSize){
    float* output = malloc(sizeof(float)*imageSize);
    
    for (int i = 0; i < imageSize; i++) {
        output[i] = input[i]*255.0f;
    }
    return output;
}

float* upcastToFloat(uint8* input, uint32 imageSize){
    float* output = malloc(sizeof(float)*imageSize);
    for (int i = 0; i < imageSize; i++) {
        output[i] = (float)input[i];
    }
    return output;
}

float* upcastToFloatAndNormalize(uint8* input, uint32 imageSize){
    float* output = malloc(sizeof(float)*imageSize);
    for (int i = 0; i < imageSize; i++) {
        output[i] = ((float)input[i])/255.0f;
    }
    return output;
}

uint8* downcastToByteAndDenormalize(float* input, uint32 imageSize){
    uint8* output = malloc(sizeof(uint8)*imageSize);
    for (int i = 0; i < imageSize; i++) {
        output[i] = input[i]*255.0f;
    }
    return output;
}

uint8* downCastToByte(float* input, uint32 imageSize){
    uint8* output = malloc(sizeof(uint8)*imageSize);
    for (int i = 0; i < imageSize; i++) {
        output[i] = (uint8)input[i];
    }
    return output;
}

void imageStatistics(uint8 * input, uint32 imageSize){
    float mean = 0;
    
    for (int i = 0; i < imageSize; i++) {
        mean += input[i];
    }
    
    printf("Image has mean value %f\n", mean/imageSize);
}

void printImageSpecs(void){
    printf("color type = %d\n", config);
    printf("bits per sample = %d\n", bitsPerSample);
    printf("number of passes = %d\n", numberOfPasses);
}

void printImage(uint8 * input, uint32 imageSize){    
    for (int i = 0; i < imageSize; i++) {
        printf("Value at %i: %hhu\n", i ,input[i]);
    }
}


//Example Usage:
//uint8* tmpVals = new uint8[4];
//tmpVals[0] = 1;
//tmpVals[1] = 2;
//tmpVals[2] = 3;
//tmpVals[3] = 4;
//
//
//float* result = multiplexToFloat(tmpVals, 1);
//
//uint8* moreTmpVals = new uint8[4];
//moreTmpVals = demultToBytes(result, 4);
//printf("float value is : %i\n", moreTmpVals[0]);

float* multiplexToFloat(uint8* data, int imageSize){
    
    float* resultingValues = malloc(sizeof(float)*imageSize);
    
    int j = 0;
    for (int i = 0; i < imageSize; i+=4) {
        union FloatAndByte aFloatAndByte;
        
        aFloatAndByte.c[0] = data[i+0];
        aFloatAndByte.c[1] = data[i+1];
        aFloatAndByte.c[2] = data[i+2];
        aFloatAndByte.c[3] = data[i+3];
        
        resultingValues[j] = aFloatAndByte.f;
        j++;
    }
    return resultingValues;
}

uint8* demultToBytes(float* data, int imageSize){
    
    uint8* resultingValues = malloc(sizeof(uint8)*imageSize);
    
    int j = 0;
    for (int i = 0; i < imageSize; i+=4) {
        union FloatAndByte aFloatAndByte;
        aFloatAndByte.f = data[j];
        
        resultingValues[i+0] = aFloatAndByte.c[0];
        resultingValues[i+1] = aFloatAndByte.c[1];
        resultingValues[i+2] = aFloatAndByte.c[2];
        resultingValues[i+3] = aFloatAndByte.c[3];
        
        j++;
    }
    return resultingValues;
}

void initializeRedTileRowPtr(void){
    rowPointers = (png_bytep*) malloc(sizeof(png_bytep) * imageLength);
    for (y=0; y<imageLength; y++){
        rowPointers[y] = (png_byte*) malloc(imageWidth*4);
    }
}

uint8* createBlackTile(void){

    uint8* image = malloc(sizeof(uint8)*10*10*4);
    imageWidth = 10;
    imageLength = 10;
    
    initializeRedTileRowPtr();

    config = 6;
    bitsPerSample = 8; 
    numberOfPasses = 1;
    
    imageWidth = imageWidth;
    imageLength = imageLength;

    bitsPerSample = bitsPerSample; // = 8 bits
    samplesPerPixel = 4; // = 4 bytes
    
    bitsPerPixel = samplesPerPixel*bitsPerSample;
    linebytes = samplesPerPixel * imageWidth; // = 640
    
    bitsPerPixel = samplesPerPixel*bitsPerSample;
    linebytes = samplesPerPixel * imageWidth; // = 640
    //linebytes = png_get_rowbytes(pngPtr, infoPtr); = 640
    imageBitSize = (sizeof(uint8) * imageWidth * imageLength * samplesPerPixel);
    imageSize = imageWidth * imageLength * samplesPerPixel;
    //printf("linebytes = %i, expected %i\n",linebytes,png_get_ro
    
    
    for (int y = 0; y<10; y++) {
        for (int x = 0; x<40; x+=4) {
            image[(y*40)+x+0]=0;
            image[(y*40)+x+1]=0;
            image[(y*40)+x+2]=0;
            image[(y*40)+x+3]=255;
        }
    }
    return image;
}

uint32 getImageSize(void){
    return imageSize;
}

uint32 getImageSizeInFloats(void){
    return imageSize*sizeof(float);
}

uint32 getImageLength(void){
    return imageLength;
};

uint32 getImageHeight(void){
    return imageLength;
};

uint32 getImageWidth(void){
    return imageWidth;
};

uint32 getConfig(void){
    return config;
};

uint32 getBitsPerSample(void){
    return bitsPerSample;
};

uint32 getSamplesPerPixel(void){
    return samplesPerPixel;
};

uint32 getImageRowPitch(void){
    return samplesPerPixel * imageWidth;
};

