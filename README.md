# A Two level Spatial In Memory Index for Non point-Data

Source code for the 'A Two-level Spatial In-Memory Index for Non-point Data' paper submitted to SIGMOD 2021 with paper ID ##.

### Dependencies
- g++/gcc
- Boost Library 
- OpenMP

**Note** that Boost library is required for all the implementations using R-Tree.

### Datasets

A sample of the datasets used can be found [here](https://1drv.ms/u/s!AkJtflSedw1rdW0uoOSJt_zOF9k?e=PrZE70). 
Files that include ```_mbr``` contain only the MBRs of the objects, while files that include ```_geom``` contain the full geometry of the objects and should be used when running refinement experiments. Files that include ```_90%``` and ```_10%``` should be used when running update experiments. Finally, there are three query sets. ```TIGER_c0.1%_n10000.qry``` should be used with TIGER dataset. Equally, ```uni_c0.1%_n10000.qry```  and ```zipf_c0.1%_n10000.qry``` should be used with uniform and zipfian datasets respectively.

**Important:** When using Quad-Tree with zipfian dataset you must limit the height of the tree using ```-h```.


### Compile
Compile using ```make all``` or ```make <option>``` where <option> can be one of the following:
   - oneLevel 
   - twoLevel 
   - twoLevelPlus 
   - rtree 
   - quadTree 
   - refinement 
   - batch 
   - update

### Parameters
We tried to be as consistent as possible when it comes to parameters. As a result, all the programs use the following parameters:
| Parameters | README |
| ------ | ------ |
| -p | The number of partitions to be used |
| -w | Window query |
| -d | Disk Query. **Keep in mind** that when using -d the use of -e is mandatory. (see examples) |
| -e | Radius of the disk query |
| -v | Specify different versions. Used only in refinement.cpp and batchProcessing.cpp. (see examples) |
| -i | Number of iterations |
| -c | Set the capacity of a quadrant. Only used with QuadTree|
| -h |  Set the height of the QuadTree. Only used with QuadTree|
| -m |  Print the manual|

### Files

- ##### partition/partition.h: 
    
    Contains all the required partition methods for the different implemenations used in our experiments.  

- ##### main_onelevel.cpp - grid/oneLevel.h: 
    
    Contains the code used for the 1-level experiments. 1-level supports both window and disk queries.
    ##### Example 
    - ###### window Query
    
        ```sh
        $ ./oneLevel -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### disk Query
    
        ```sh
        $ ./oneLevel -p 3000 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_twolevel.cpp - grid/twoLevel.h: 
    
    Contains the code used for the 2-level experiments. 2-level supports both window and disk queries.
    ##### Example 
    - ###### window Query
    
        ```sh
        $ ./twoLevel -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### disk Query
    
        ```sh
        $ ./twoLevel -p 3000 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_twolevel_plus.cpp - grid/twoLevelPlus.h: 
    
    Contains the code used for the 2-level+ experiments. 2-level+ supports only window queries.
    ##### Example 
    - ###### window Query
    
        ```sh
        $ ./twoLevelPlus -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_rtree.cpp: 
    
    Contains the code used for the R-Tree experiments. R-Tree supports both window and disk queries. Instead of partitions, R-Tree uses maximum node capacity. In our experiments we used 16, but the value can be changed by changing the value of  ``` #define MAX_NODE_CAPACITY 16``` to the desired value.
     ##### Example 
    - ###### window Query
    
        ```sh
        $ ./rtree -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### disk Query
    
        ```sh
        $ ./rtree -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_quadtree.cpp - quadTree/: 
    Contains the code used for the Quad-Tree experiments. Quad-Tree supports both window and disk queries. Parameter ```-c``` is required, while parameted ```-h``` can be omitted. Based on our experiments, we strongly suggest that the capacity is great than 100. 
    ##### Example 
    - ###### window Query with quadrant capacity of 1000
    
        ```sh
        $ ./qt -c 1000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### disk Query with quadrant capacity of 1000 and maximum height of 10
    
        ```sh
        $ ./qt -c 1000 -h 10 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
- ##### main_refine.cpp - grid/refinement.h: 
    Contains the code used for the refinement experiments. Refinement supports both window and disk queries. As described in the paper we experimented with three different approaches; *Simple*, *RefAvoid*, *RefAvoid+* (see figure 7). Window queries can perform all three methods, while disk queries can only perform *Simple* and *RefAvoid*. Using parameter ```-v``` the desired method can be specified. Valid values are:
    - 1 for *Simple*
    - 2 for *RefAvoid*
    - 3 for *RefAvoid+*
    
    ##### Example 
    - ###### window Query with *RefAvoid*
    
        ```sh
        $ ./refine -p 3000 -w -v 2 TIGER_ROADS_geom.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### disk Query with *RefAvoid*
    
        ```sh
        $ ./refine -p 3000 -d -e 0.1 -v 2 TIGER_ROADS_geom.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_batch.cpp - grid/batchProcessing.h: 
    Contains the code used for the batch processing experiments. Batch processing supports only window queries with the two different approaches described in the paper; *Queries-Based* and *Tiles-based*. Using parameter ```-t``` the number of threads can be specified. Again using parameter ```-v``` the desired approach can be specified. Valid values are:
    - 1 for *Queries-Based*
    - 2 for *Tiles-based*

    ##### Example 
    - ###### *Queries-Based* window query using 1 thread
    
        ```sh
        $ ./batch -p 3000 -t 1 -v 1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```
    - ###### *Tiles-Based* window query using 4 threads
    
        ```sh
        $ ./batch -p 3000 -t 4 -v 2 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry
        ```

- ##### main_update.cpp: 
    Contains the code used for the updates experiments. Updates can be performed when using *1-Level*, *2-Level*, *QuadTree* and *R-Tree*. Valid parameters are:
    - ```-1``` for *1-Level*
    - ```-2``` for *2-Level*
    - ```-q``` for *QuadTree*
    - ```-r``` for *R-Tree*.

    ##### Example 
    - ###### Updates on QuadTree with capacity 1000
    
        ```sh
        $ ./update -c 1000 -q TIGER_ROADS_mbr_90%.inp TIGER_c0.1%_n10000.qry TIGER_ROADS_mbr_10%.inp
        ```
    - ###### Update on *2-Level*
    
        ```sh
        $ ./update -p 3000 -2 TIGER_ROADS_mbr_90%.inp TIGER_c0.1%_n10000.qry TIGER_ROADS_mbr_10%.inp
        ```
