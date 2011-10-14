/*
 *  FileHandler.h
 *  MedialAxisTransform
 *
 *
 *  Created by Beau Johnston on 25/08/11.
 *  Copyright (C) 2011 by Beau Johnston.
 *
 *  Please email me if you have any comments, suggestions or advice:
 *                              beau@inbeta.org
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


#ifndef openCLImageLoad_FileHandler_h
#define openCLImageLoad_FileHandler_h

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DymanicArray.h"

    
#define bool int 
#define true 1 
#define false 0
//#define NULL ((void *)0) 
    
    /* private functions */
    int getdir (char* dir, char* ext);
    int getFilesInDirectoryWithName (char* dir, char* name);
    char* removeSubstring(char* str, char* substr);
    bool doesNumberingStartAt0(char* name, char* ext);
    bool doesNumberingStartAt1(char* name, char* ext);
    bool noFilesAreMissing(char* name, char* ext);

    size_t locateLast(char* string, char pattern);
    char *replace_str(const char *str, const char *orig, const char *rep);
    char* substring(size_t i, size_t j, char *ch);
    void generateListOfAssociatedFiles(char* filename);
    char* getNextFileName(void);
    bool areFilesLeft(void);
    int numberOfFiles(void);
    void printFiles(void);
    void sortFilesNumerically(char* name, char* ext);
    char* nullTerminateString(const char * in);
    
#endif
