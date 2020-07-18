#include <vector>
#include <iostream>
#include "./containers/relation.h"
#include "def.h"
#include "./quadTree/QuadTree.h"
#include "getopt.h"

using namespace std;

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./qt - code for the Quad-Tree processing experiments" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./qt [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       Mandatory arguments" << endl << endl;
    cerr << "       -c" << endl;
    cerr << "              node capacity" << endl;
    cerr << "       -w" << endl;
    cerr << "              window query" << endl;
    cerr << "       -d" << endl;
    cerr << "              disk query. Option -e should be used when using -d " << endl;
    cerr << "       -e" << endl;
    cerr << "              radius of the disk query" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -h" << endl;
    cerr << "              heigh of Quad-Tree" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "        Window query using Quad-Tree with capacity of 1000." << endl;
    cerr << "              ./qt -c 1000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "        Disk query using Quad-Tree with capacity of 1000, maximum height of 10 and radius of 0.1" << endl;
    cerr << "              ./qt -c 1000 -h 10 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    Relation R, S;
    unsigned long long result = 0;
    unsigned long long result1 = 0;
    Timer tim;
    char *k;
    int method = -1;
    int runNumThreads = 1;
    double timeIndexingOrPartitioning = 0, queryExecutionTime = 0, timeCounting = 0;
    char c;
    int  NODE_CAPACITY = -1;
    int level = -1;
    int queryCase, queryMethod=-1;
    Coord epsilon = -1.0;
    Coord minX, maxX, minY, maxY;
    Coord diffX, diffY, maxExtend;
    QuadTreeNode* tree;
    vector<RecordId> resultItemsIds;
    int NUM_ITERATIONS = 1;


    while ((c = getopt(argc, argv, "c:h:e:wdi:m")) != -1)
    {
        switch (c)
        {
            case 'c':
                NODE_CAPACITY = atoi(optarg);
                break;          
            case 'h':
                level = atoi(optarg);
                break;
            case 'w':
                queryMethod = WINDOW_QUERY;
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'd':
                queryMethod = DISK_QUERY;
                break;
            case 'i':
                NUM_ITERATIONS = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cout << "Wrong arguments! ";
                exit(0);
                break;
        }
    }

    if(NODE_CAPACITY == -1)
    {
        cerr << "Node capacity is missing" << endl;
        usage();
    }

    if(queryMethod == -1)
    {
        cerr << "Query method is missing" << endl;
        usage();
    }

    if (queryMethod == DISK_QUERY){
        if (epsilon < 0){
            cout<<"Radius value is missing"<<endl;
            usage();
        }
    }
    
    if (NODE_CAPACITY < 100){
        cout<<"Capacity value should be greater than 100"<<endl;
        usage();
    }

    //Load inputs
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            R.load(argv[optind]);
        }
        #pragma omp section
        {
            S.load(argv[optind+1]);
        }
    }

    BoundingBox boundary(0.0, 0.0,1.0,1.0);

    cout << "Method\t" << "Query ID\t" << "Total Results\t" << "Indexing Time\t" << "Total Filtering Time" << endl; 


    switch (queryMethod){
        case WINDOW_QUERY: 
            minX = min(R.minX, S.minX);
            maxX = max(R.maxX, S.maxX);
            minY = min(R.minY, S.minY);
            maxY = max(R.maxY, S.maxY);
            diffX = maxX - minX;
            diffY = maxY - minY;
            maxExtend = (diffX<diffY)?diffY:diffX;

            R.normalize(minX, maxX, minY, maxY, maxExtend);
            S.normalize(minX, maxX, minY, maxY, maxExtend);
            if( level == -1 ){
                tim.start();
                tree = new QuadTreeNode(boundary, NODE_CAPACITY, 0);
                tree->loadDataNoHeight(R);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) {
                        auto &s = S[i];

                        tree->rangeQuery(s,resultItemsIds);

                        result1 = resultItemsIds.size();
                        resultItemsIds.clear();
                        resultItemsIds.reserve(result1);
                        
                        tim.start();
                        tree->rangeQuery(s,resultItemsIds);
                        queryExecutionTime = tim.stop();

                        unsigned long long result = 0;
                        tim.start();
                        for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                            result++;
                        }
                        timeCounting = tim.stop();
                        resultItemsIds.clear();
                        cout << "QuadTreeNoHeight-Window-Query\t"<<i << "\t" << result <<"\t" << timeIndexingOrPartitioning <<"\t"<<queryExecutionTime+timeCounting<< endl;
                    }
                }

            } 
            else
            {
                tim.start();
                tree = new QuadTreeNode(boundary, NODE_CAPACITY, 0, level);
                tree->loadData(R);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) {
                        auto &s = S[i];

                        tree->rangeQuery(s, resultItemsIds);

                        result1 = resultItemsIds.size();
                        resultItemsIds.clear();
                        resultItemsIds.reserve(result1);
                        
                        tim.start();
                        tree->rangeQuery(s,resultItemsIds);
                        queryExecutionTime = tim.stop();

                        unsigned long long result = 0;
                        tim.start();
                        for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                            result++;
                        }
                        timeCounting = tim.stop();

                        resultItemsIds.clear();
                        cout << "QuadTree-Window-Query\t"<<i<<"\t" << result <<"\t" << timeIndexingOrPartitioning <<"\t"<<queryExecutionTime+timeCounting<< endl;
                    }
                }
            }
        break;

        case DISK_QUERY:
            minX = R.minX;
            maxX = R.maxX;
            minY = R.minY;
            maxY = R.maxY;
            diffX = maxX - minX;
            diffY = maxY - minY;
            maxExtend = (diffX<diffY)?diffY:diffX;
            
            R.normalize(minX, maxX, minY, maxY, maxExtend);
            S.normalize(minX, maxX, minY, maxY, maxExtend);
            
            double xMin = -117;
            double yMin = 33;
            double xMax = -78;
            double yMax = 46;

            xMin = Coord(xMin - minX) / maxExtend;
            xMax   = Coord(xMax   - minX) / maxExtend;
            yMin = Coord(yMin - minY) / maxExtend;
            yMax   = Coord(yMax   - minY) / maxExtend;

            Coord world = (xMax - xMin) * (yMax - yMin); 
            
            epsilon = sqrt(epsilon*world/3.14);

            S.loadDisk(epsilon);

            if( level == -1 ){
                tim.start();
                tree = new QuadTreeNode(boundary, NODE_CAPACITY, 0);
                tree->loadDataNoHeight(R);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) {
                        auto &s = S[i];

                        tree->rangeQuery(s, resultItemsIds, epsilon);

                        result1 = resultItemsIds.size();
                        resultItemsIds.clear();
                        resultItemsIds.reserve(result1);
                        
                        tim.start();
                        tree->rangeQuery(s,resultItemsIds,epsilon);
                        queryExecutionTime = tim.stop();

                        unsigned long long result = 0;
                        tim.start();
                        for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                            result++;
                        }
                        timeCounting = tim.stop();

                        resultItemsIds.clear();

                        cout << "QuadTreeNoHeight-Disk-Query\t"<<i<<"\t" << result <<"\t" << timeIndexingOrPartitioning <<"\t"<<queryExecutionTime+timeCounting<< endl;
                    }
                }

            } 
            else
            {
                tim.start();
                tree = new QuadTreeNode(boundary, NODE_CAPACITY, 0, level);
                tree->loadData(R);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) {
                        auto &s = S[i];

                        tree->rangeQuery(s, resultItemsIds, epsilon);

                        result1 = resultItemsIds.size();
                        resultItemsIds.clear();
                        resultItemsIds.reserve(result1);
                        
                        tim.start();
                        tree->rangeQuery(s,resultItemsIds,epsilon);
                        queryExecutionTime = tim.stop();

                        unsigned long long result = 0;
                        tim.start();
                        for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                            result++;
                        }
                        timeCounting = tim.stop();

                        resultItemsIds.clear();
                        cout << "QuadTree-Disk-Query\t"<<i<< "\t" << result <<"\t" << timeIndexingOrPartitioning <<"\t"<<queryExecutionTime+timeCounting<< endl;
                    }
                }
            }

        break;
    }
    
    tree->deleteTree();

    return 0;
}

