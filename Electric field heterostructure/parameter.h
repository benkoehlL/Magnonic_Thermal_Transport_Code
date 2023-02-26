#ifndef PARAMETER_H
#define PARAMETER_H

#define J 1
#define S 0.5
#define Dmax 1
#define Dincrement 0.00001
#define freefieldsize 1
#define fieldsize 46//Due to periodic boundaries use even values only
#define xsize (2*freefieldsize+fieldsize)
#define ysize 30//Due to periodic boundaries use even values only
#define systemsize xsize*ysize
#define epsilon 0.0000001
#define cutoffdelta 0.005
#define cutoffdeltados 0.01
#endif
