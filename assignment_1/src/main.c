/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "affinity.h"
#include "allocate.h"
#include "timing.h"

int main(int argc, char** argv) {
	const char* node_env = getenv("SLURM_JOB_NODELIST");
	const char* name = "World";
	if(argc>=1){
		name = argv[1];
	}
	printf("Hello , %s , ! \n",name);       
	if(node_env){
		printf("The node : %s\n",node_env);
	}
	return EXIT_SUCCESS; }
