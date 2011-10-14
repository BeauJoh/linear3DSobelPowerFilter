/*
 *  FileHandler.c
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


#include "FileHandler.h"

/* private variables */
DynamicArray files;
DynamicArray orderedFiles;
//char* files[1000];
int fileSize=0;
//char * orderedFiles[1000];
int traverser;
char* path;



int getdir (char* dir, char* ext){
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir)) == NULL) {
        printf("error (%i) opening %s \n",errno,dir);        
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (strstr(dirp->d_name,ext) != NULL) {
            addElement(&files, dirp->d_name);
            //files[fileSize] = dirp->d_name;
            fileSize++;
        }
    }
    closedir(dp);
    return 0;
}

int getFilesInDirectoryWithName (char* dir, char* name){
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir)) == NULL) {
        printf("error (%i) opening %s \n",errno,dir);        
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (strstr(dirp->d_name, name) != NULL) {
            //files[fileSize]= malloc(sizeof(char)*1000);
            //memset(files[fileSize], '\0', 1000);
            
            //strcpy(files[fileSize],dirp->d_name);
            fileSize++;
            addElement(&files, dirp->d_name);
        }
    }
    closedir(dp);
    return 0;
}

//void escape_str(char *dest, char *src)
//{
//    *dest = 0;
//    while(*src)
//    {
//        switch(*src)
//        {
//            case '\n' : strcat(dest++, "\\n"); break;
//            case '\"' : strcat(dest++, "\\\""); break;
//            default:  *dest = *src;
//        }
//        *src++;
//        *dest++;
//        *dest = 0;                     
//    }
//}

char * replace_str(
               char const * const original, 
               char const * const pattern, 
               char const * const replacement
               ) {
    size_t const replen = strlen(replacement);
    size_t const patlen = strlen(pattern);
    size_t const orilen = strlen(original);
    
    size_t patcnt = 0;
    const char * oriptr;
    const char * patloc;
    
    // find how many times the pattern occurs in the original string
    for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
    {
        patcnt++;
    }
    
    {
        // allocate memory for the new string
        size_t const retlen = orilen + patcnt * (replen - patlen);
        char * const returned = (char *) malloc( sizeof(char) * (retlen + 1) );
        
        if (returned != NULL)
        {
            // copy the original string, 
            // replacing all the instances of the pattern
            char * retptr = returned;
            for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
            {
                size_t const skplen = patloc - oriptr;
                // copy the section until the occurence of the pattern
                strncpy(retptr, oriptr, skplen);
                retptr += skplen;
                // copy the replacement 
                strncpy(retptr, replacement, replen);
                retptr += replen;
            }
            // copy the rest of the string.
            strcpy(retptr, oriptr);
        }
        return returned;
    }
}


//char *replace(char *str, char *orig, char *rep)
//{
//    char* buffer = malloc(sizeof(char)*4096);
//    //memset(buffer, 0, 4096);
//
//    char *p;
//    
//    if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
//        return str;
//        
//    strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
//    //buffer[p-str] = '\0';
//    //escape_str(buffer, buffer);
//
//    sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));
//    return buffer;
//}

char* nullTerminateString(const char * in){
    printf("strlen(in) %lu\n", strlen(in));
    char* resultString = malloc(sizeof(char*)*(strlen(in)+1));
    strcpy(resultString, in);
    strcat(resultString, "\0");
    return resultString;
}

char* removeSubstring(char* str, char* substr){
    return replace_str(str, substr, "");
}

bool doesNumberingStartAt0(char* name, char* ext){
    int lowest = numberOfFiles();
    
    for (int j=0; j < numberOfFiles(); j++) {
        int i;
        
        char* tmp = removeSubstring(getElementAtIndex(files, j), name);
        char* inNumber = removeSubstring(tmp, ext);
        
        i = atoi(inNumber);
        if(i)
        {
            if (i < lowest) {
                lowest = i;
            }
        }
    }
    
    if (lowest == 0) {
        return true;
    }
    else{
        return false;
    }

}


bool doesNumberingStartAt1(char* name, char* ext){
    int lowest = numberOfFiles();
    
    for (int j=0; j < numberOfFiles(); j++) {
        int i;
        char* tmp = removeSubstring(getElementAtIndex(files, j), name);
        char* inNumber = removeSubstring(tmp, ext);
        
        i = atoi(inNumber);
        if(i)
        {
            if (i < lowest) {
                lowest = i;
            }
        }
    }
    
    if (lowest == 1) {
        return true;
    }
    else{
        return false;
    }
}


bool noFilesAreMissing(char* name, char* ext){
    int lowest = numberOfFiles();
    int highest = 0;
    int tempVec[numberOfFiles()];
    int tempVecSize = 0;
    for (int j=0; j < numberOfFiles(); j++) {
        int i;
        char* tmp = removeSubstring(getElementAtIndex(files, j), name);
        char* inNumber = removeSubstring(tmp, ext);
        
        i = atoi(inNumber);
        if(i)
        {
            if (i < lowest) {
                lowest = i;
            }
            if (i > highest) {
                highest = i;
            }
            tempVec[j] = i;
            tempVecSize ++;
        }
    }
    
    int tempArray[highest-lowest];
    //populate with garbage
    for (int j = 0; j < tempVecSize; j++) {
        tempArray[j] = -1;
    }
    
    //populate with real values
    for (int j = 0; j < tempVecSize; j++) {
        tempArray[tempVec[j]] = tempVec[j];
    }
    
    for (int k = lowest; k < highest; k++) {
        if(tempArray[k]==-1){
            printf(" Image number : %i dosn't exist! \n",k);
            printf(" It is compulsory that the image stack be fully populated\n");
            printf(" Please fix this and try again!\n");
            return false;
        }
    }
    return true;    
}

char* substring(size_t i, size_t j, char *ch) {
    if (i > j) {
        printf("starting index cannot be greated than finishing\n");
        return '\0';
    }
    //if (i < 0 || j < 0) {
    //    printf("Negative indexing forbidden\n");
    //    return '\0';
    //}
    
    char* resultString = malloc(sizeof(char)*(j-i));
    
    for (int k = 0; k < j-i; k++) {
        resultString[k] = ch[i + k];
    }
    
    //null terminate string
    resultString[j-i]='\0';
    
    return resultString;
}

size_t locateLast(char* string, char pattern){
    return strrchr(string, pattern)-string;
}

void generateListOfAssociatedFiles(char* filename){

    //initialise dynamic arrays
    files = initDynamicArray();
    orderedFiles = initDynamicArray();
    
    char* handedInString = filename;
    
    char* file = substring(locateLast(handedInString,'/')+1, strlen(handedInString), handedInString);
    
    
    path = substring(0, locateLast(handedInString,'/')+1, handedInString);

    char* cutDownFile = substring(0, locateLast(file,'.'), file);
    char* extension = substring(locateLast(file,'.'), (int)strlen(file),file);
    
    printf("file has name: %s\n", file);
    printf("path has name: %s\n", path);
    printf("cut down file has name: %s\n", cutDownFile);
    printf("extension has name: %s\n", extension);
    
    char* dir = path;
    
    if (strcmp(cutDownFile, "") == 0) {
        getdir(dir, extension);
    } else {
        getFilesInDirectoryWithName(dir,cutDownFile);
    }
    
    
    traverser = 0;
    sortFilesNumerically(cutDownFile, extension);
    for (unsigned int i = 0;i < fileSize;i++) {
        printf("%s\n", getElementAtIndex(orderedFiles, i));
    }
    return;    
}

char* getNextFileName(void){
    
    char* result = malloc(sizeof(char)*strlen(getElementAtIndex(orderedFiles, traverser)
));
    strcpy(result, path);
    strcat(result, getElementAtIndex(orderedFiles, traverser)); 
    //printf("getNextFileName name result is: %s\n", getElementAtIndex(orderedFiles, traverser));
    traverser ++;
    return result;
}

bool areFilesLeft(void){
    if (traverser < fileSize+1) {
        return true;
    }
    return false; 
}

int numberOfFiles(void){
    return fileSize;
}

void printFiles(void){
    for (unsigned int i = 0; i < fileSize; i++) {
        printf("%s", getElementAtIndex(orderedFiles, i));
    }
    return;
    
}

void sortFilesNumerically(char* name, char* ext){
    if (!doesNumberingStartAt1(name, ext)&&!doesNumberingStartAt0(name, ext)){
        printf("Error, numbering for image stack doesn't start at 0 or 1\n");
        return;
    }
    
    bool startingAtZero = doesNumberingStartAt0(name, ext);
    
    
    if (!noFilesAreMissing(name, ext)) {
        printf("ERROR: files are missing, returning!\n");
        return;
    }
    
//    bool startingAtZero = false;
    
    for (int j=0; j < numberOfFiles(); j++) {
        int i;
        
        char* tmp = removeSubstring(getElementAtIndex(files, j), name);
        char* inNumber = removeSubstring(tmp, ext);
        i = atoi(inNumber);
        printf("i has value %i \n", i);

        if(i)
        {
//            cout << i << endl;
            if (startingAtZero) {
                setElementAtIndex(&orderedFiles, i, numberOfFiles(), nullTerminateString(getElementAtIndex(files, j)));
//                setElementAtIndex(&orderedFiles, i, numberOfFiles(), getElementAtIndex(files, j));
            } else{
                setElementAtIndex(&orderedFiles, i-1, numberOfFiles(), nullTerminateString(getElementAtIndex(files, j)));
//                setElementAtIndex(&orderedFiles, i-1, numberOfFiles(), getElementAtIndex(files, j));
            }
        }
    }
    //orderedFiles[numberOfFiles()] = '\0';
    
}

