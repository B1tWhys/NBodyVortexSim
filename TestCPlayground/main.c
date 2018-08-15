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

void testF(int *arrPtr[]) {
	for (int i = 0; i < 5; i++) {
		(*arrPtr)[i] = i;
	}
}

int main() {
	int *test = malloc(5 * sizeof(int));
	testF(&test);
	printf("bp\n");
}














