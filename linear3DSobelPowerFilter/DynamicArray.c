//
//  DynamicArray.c
//  fastFourierTransform
//
//  Created by Beau Johnston on 26/09/11.
//  Copyright 2011 University Of New England. All rights reserved.
//

#include "DymanicArray.h"

DynamicArray initDynamicArray(void){   
    DynamicArray dynamicArray;
    
    dynamicArray.theArray = NULL;
    dynamicArray.numElements = 0; //number of elements being used
    dynamicArray.numAllocated = 0; //size of the array allocated
    return dynamicArray;
}

void destroyDynamicArray(DynamicArray * dynamicArray){
    free ((*dynamicArray).theArray);
}

void setElementAtIndex(DynamicArray * dynamicArray, int i, int sizeOfArray, CString item){
    //invalid parameters if the index to add is larger than the expected array size
    assert(i < sizeOfArray);
    
    //array already fully populated
    assert((*dynamicArray).numElements < sizeOfArray);
    
    //memory allocation stuff
    if ((*dynamicArray).numAllocated < sizeOfArray) {
        (*dynamicArray).theArray = malloc(sizeof(CString)*sizeOfArray);
        (*dynamicArray).numAllocated = sizeOfArray;
    }
    
    printf("item has value: %s\n", item);
    
    //add element stuff
    (*dynamicArray).theArray[i] = item;
    (*dynamicArray).numElements ++;
}

void addElement(DynamicArray * dynamicArray, CString item){
    if((*dynamicArray).numElements == (*dynamicArray).numAllocated) // Are more refs required?
    {
        if ((*dynamicArray).numAllocated == 0)
            (*dynamicArray).numAllocated = 1; // Start off with 1 refs
        else
            (*dynamicArray).numAllocated *= 2; // Double the number
        // of refs allocated
        
        // Make the reallocation transactional
        // by using a temporary variable first
        void *_tmp = realloc((*dynamicArray).theArray, ((*dynamicArray).numAllocated * sizeof(CString)));
        
        // If the reallocation didn't go so well,
        // inform the user and bail out
        if (!_tmp)
        {
            //printf(stderr, "ERROR: Couldn't realloc memory!\n");
            return;
        }
        
        // Things are looking good so far
        (*dynamicArray).theArray = (CString*)_tmp;
    }
    
    (*dynamicArray).theArray[(*dynamicArray).numElements] = item;
    (*dynamicArray).numElements++;

    return;
}

CString getElementAtIndex(DynamicArray dynamicArray, int i){
    //printf("%s", dynamicArray.theArray[i]);
    return dynamicArray.theArray[i];
}

CString getLastElementAndRemove(DynamicArray * dynamicArray){
    (*dynamicArray).numElements --;
    return (*dynamicArray).theArray[(*dynamicArray).numElements+1];
}


int getNumberOfElements(DynamicArray dynamicArray){
    return dynamicArray.numElements;
}

