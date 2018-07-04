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
	int mod1 = -4;
	int mod2 = 10;
	int res = (mod1 < 0) ? mod2 + (mod1 % mod2) : mod1 % mod2;
	
	printf("%i\n", res);
}

#pragma mark - main

int main(int argc, const char * argv[]) {
	testFunc();
	
	return 0;
}
















