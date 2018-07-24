//
//  main.c
//  TestCPlayground
//
//  Created by Skyler Arnold on 6/13/18.
//  Copyright © 2018 Skyler Arnold. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

void testFunc() {
	FILE *testFile;
	testFile = fopen("./testFile", "w");
	
	
//	fprintf(testFile, "test text 123");
	fputc(29, testFile);
	
}

#pragma mark - main

int main(int argc, const char * argv[]) {
	testFunc();
	
	return 0;
}
















