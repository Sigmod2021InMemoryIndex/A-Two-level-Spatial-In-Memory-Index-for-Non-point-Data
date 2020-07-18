#pragma once
#ifndef _DEF_H_
#define _DEF_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "omp.h"

using namespace std;

#define WINDOW_QUERY                                                    0
#define DISK_QUERY                                                      1

#define QUERIES_BASED                                                   1
#define TILES_BASED                                                     2

#define SIMPLE_REFINEMENT                                               1
#define REFAVOID                                                        2
#define REFAVOIDPLUS                                                    3

#define R_TREE                                                          1
#define QUAD_TREE                                                       2
#define ONE_LEVEL                                                       3
#define SECOND_LEVEL                                                    4

#define EPS 1e-08

typedef double Coord;
typedef size_t RecordId;

class Record;
class Relation;


class Timer{
    private:
        using Clock = std::chrono::high_resolution_clock;
        Clock::time_point start_time, stop_time;

    public:
        Timer()
        {
                start();
        }

        void start()
        {
                start_time = Clock::now();
        }

        double getElapsedTimeInSeconds()
        {
                return std::chrono::duration<double>(stop_time - start_time).count();
        }

        double stop()
        {
                stop_time = Clock::now();
                return getElapsedTimeInSeconds();
        }
    };


    
    

#endif
