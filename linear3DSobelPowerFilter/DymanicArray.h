//
//  DymanicArray.h
//  fastFourierTransform
//
//  Created by Beau Johnston on 26/09/11.
//  Copyright 2011 University Of New England. All rights reserved.
//

#ifndef fastFourierTransform_DymanicArray_h
#define fastFourierTransform_DymanicArray_h

#include <stdlib.h>
#include <assert.h>
#include <string.h>


typedef char * CString;

typedef struct __DynamicArray{
    CString *theArray;
    int numElements; //number of elements being used
    int numAllocated; //size of the array allocated
}DynamicArray;

DynamicArray initDynamicArray(void);
void destroyDynamicArray(DynamicArray * dynamicArray);

CString getElementAtIndex(DynamicArray dynamicArray, int i);
void setElementAtIndex(DynamicArray * dynamicArray, int i, int sizeOfArray, CString item);

void addElement(DynamicArray * dynamicArray, CString item);
CString getElement(DynamicArray * dynamicArray);

int getNumberOfElements(DynamicArray dynamicArray);
CString getLastElementAndRemove(DynamicArray * dynamicArray);


#endif
