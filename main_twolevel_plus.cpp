#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./partitioning/partition.h"
#include "./grid/twoLevelPlus.h"

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./twoLevelPlus - range query using the 2-level+ algorithm" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./twoLevelPlus [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       Mandatory arguments" << endl << endl;
    cerr << "       -p" << endl;
    cerr << "              number of partions per dimension" << endl;
    cerr << "       -w" << endl;
    cerr << "              window query" << endl;
    cerr << "       -d" << endl;
    cerr << "              disk query. Option -e should be used when using -d " << endl;
    cerr << "       -e" << endl;
    cerr << "              radius of the disk query" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "       Window range query using the 2-level+ algorithm with 3000 partitions per dimension." << endl;
    cerr << "              ./twoLevelPlus -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    char c;
    int runNumPartitionsPerRelation = -1;
    Timer tim;
    double timeCounting = 0;
    double timeSorting = 0, timeIndexingOrPartitioning = 0, timeFindingTiles = 0, timeCorner = 0,timeInside = 0,timeBorders = 0;
    Relation R, S, *pR;
    vector<Decompose> *pRXStart , *pRXEnd, *pRYStart, *pRYEnd;
    size_t *pRA_size, *pRB_size, *pRC_size, *pRD_size;
    int runNumPartitions = -1;
    vector<Coord> indexR;
    vector<Coord> indexS;
    int queryCase, queryMethod=-1;
    double epsilon = 0.0;
    vector<RecordId> resultItemsIds;
    int NUM_ITERATIONS = 1;

    unsigned long long result = 0;

    while ((c = getopt(argc, argv, "p:wdue:i:m")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
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
                NUM_ITERATIONS = atof(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cerr << "Wrong arguments! ";
                exit(1);
                break;
        }
    }

    if(runNumPartitionsPerRelation == -1)
    {
        cerr << "Number of partitions is missing" << endl;
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

    // Load inputs
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

    Coord minX = min(R.minX, S.minX);
    Coord maxX = max(R.maxX, S.maxX);
    Coord minY = min(R.minY, S.minY);
    Coord maxY = max(R.maxY, S.maxY);
    Coord diffX = maxX - minX;
    Coord diffY = maxY - minY;
    Coord maxExtend = (diffX<diffY)?diffY:diffX;

    R.normalize(minX, maxX, minY, maxY, maxExtend);
    S.normalize(minX, maxX, minY, maxY, maxExtend);

    indexR.resize(4*R.numRecords);
    indexS.resize(4*S.numRecords);
    
    int counter = 0;
    for ( int i = 0 ; i < R.numRecords ; i ++ ){
        
        indexR[counter] = R[i].xStart;
        counter++;
        indexR[counter] = R[i].xEnd;
        counter++;
        indexR[counter] = R[i].yStart;
        counter++;
        indexR[counter] = R[i].yEnd;
        counter++;
    }
    counter = 0;
    for ( int i = 0 ; i < S.numRecords ; i ++ ){
        indexS[counter] = S[i].xStart;
        counter++;
        indexS[counter] = S[i].xEnd;
        counter++;
        indexS[counter] = S[i].yStart;
        counter++;
        indexS[counter] = S[i].yEnd;
        counter++;
    }

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;

    cout << "Method\t" << "Query ID\t" << "Results\t" << "Indexing Time\t" << "Filtering Time" << endl;


    
    switch (queryMethod)
    {
        case WINDOW_QUERY:
            pR = new Relation[runNumPartitions];

            pRXStart = new vector<Decompose>[runNumPartitions];
            pRXEnd = new vector<Decompose>[runNumPartitions];
            pRYStart = new vector<Decompose>[runNumPartitions];
            pRYEnd = new vector<Decompose>[runNumPartitions];

            pRA_size = new size_t[runNumPartitions];
            pRB_size = new size_t[runNumPartitions];
            pRC_size = new size_t[runNumPartitions];
            pRD_size = new size_t[runNumPartitions];
            memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRD_size, 0, runNumPartitions*sizeof(size_t));
            
            tim.start();
            partition::twoLevelPlus::PartitionTwoDimensional(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, pRXStart, pRXEnd,  pRYStart, pRYEnd, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();

            // tim.start();
            twoLevelPlus::window::sort::SortXStart(pRXStart, pRA_size , runNumPartitions);
            twoLevelPlus::window::sort::SortXEnd(pRXEnd, pRA_size , pRB_size, pRC_size, runNumPartitions);
            twoLevelPlus::window::sort::SortYStart(pRYStart, pRA_size , runNumPartitions);
            twoLevelPlus::window::sort::SortYEnd(pRYEnd, pRA_size ,pRB_size, runNumPartitions);
            // timeSorting = tim.stop();
            timeIndexingOrPartitioning = tim.stop();

            for (size_t j = 0; j < S.numRecords; j++)
            {
                for(int o = 0; o < NUM_ITERATIONS; o++) 
                {

                    auto &s = S[j];
                    result = 0;
                    vector<int> insideCells;
                    vector<int> cornerCells;

                    tim.start();
                    queryCase = twoLevelPlus::FindRelevantTiles(s, cornerCells, insideCells, runNumPartitionsPerRelation);
                    timeFindingTiles = tim.stop();

                    for(int pid = 0; pid< insideCells.size();pid++)
                    {
                        result += pRA_size[insideCells[pid]];
                    }
                    
                    for(int i = 0; i < 4;i++)
                    {
                        if(cornerCells[i] != -1)
                            twoLevelPlus::window::RangeCorners(pR[cornerCells[i]], s, 0, pRA_size[cornerCells[i]], resultItemsIds);
                    }

                    for(int i = 0; i < 2;i++)
                    {
                        if(cornerCells[i] != -1)
                            twoLevelPlus::window::RangeCorners(pR[cornerCells[i]], s, pRA_size[cornerCells[i]], pRB_size[cornerCells[i]], resultItemsIds);
                    }
                    for(int i = 0; i < 3;i+=2)
                    {
                        if(cornerCells[i] != -1)
                            twoLevelPlus::window::RangeCorners(pR[cornerCells[i]], s, pRB_size[cornerCells[i]], pRC_size[cornerCells[i]], resultItemsIds);
                    }
                    twoLevelPlus::window::RangeCorners(pR[cornerCells[0]], s, pRC_size[cornerCells[0]], pRD_size[cornerCells[0]], resultItemsIds);
                    
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevelPlus::window::RangeCorners(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], pRA_size[i], pRB_size[i],s.yStart,resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeCorners(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], pRB_size[i], pRC_size[i], s.xStart,resultItemsIds);
                            }

                            break;
                        case 3:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], 0, pRA_size[i],s.yStart,resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], pRA_size[i], pRB_size[i],s.yStart,resultItemsIds);
                            }

                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], 0, pRA_size[i],s.xStart,resultItemsIds);  
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], pRB_size[i], pRC_size[i], s.xStart,resultItemsIds);
                            }

                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch2(pRXStart[i], 0, pRA_size[i],s.xEnd,resultItemsIds);
                            }

                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch2(pRYStart[i], 0, pRA_size[i],s.yEnd,resultItemsIds);
                            }
                            break;
                    }
                    
                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);

                    tim.start();
                    auto lsize = insideCells.size();
                    for(int pid = 0; pid< lsize; pid++)
                    {
                        auto cid = insideCells[pid];
                        for (auto it = 0; it < pRA_size[cid]; it++) {
                            resultItemsIds.push_back(pRXStart[cid][it].id);
                        }
                    }
                    timeInside = tim.stop();
                    // cout << "Inside " << resultItems.size() << endl;
                    
                    tim.start();
                    for(int i = 0; i < 4;i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevelPlus::window::RangeCorners(pR[cid], s, 0, pRA_size[cid],resultItemsIds);
                    }

                    for(int i = 0; i < 2;i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevelPlus::window::RangeCorners(pR[cid], s, pRA_size[cid], pRB_size[cid],resultItemsIds);
                    }
                    for(int i = 0; i < 3;i+=2)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevelPlus::window::RangeCorners(pR[cid], s, pRB_size[cid], pRC_size[cid],resultItemsIds);
                    }
                    auto cid = cornerCells[0];
                    twoLevelPlus::window::RangeCorners(pR[cid], s, pRC_size[cid], pRD_size[cid],resultItemsIds);
                    timeCorner = tim.stop();

                    tim.start();
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevelPlus::window::RangeBClass(pRXStart[i],indexR,s, 0, pRA_size[i],resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], pRA_size[i], pRB_size[i],s.yStart,resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeCClass(pRXStart[i],indexR,s, 0, pRA_size[i],resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], pRB_size[i], pRC_size[i], s.xStart,resultItemsIds);
                            }

                            break;
                        case 3:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], 0, pRA_size[i],s.yStart,resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRYEnd[i], pRA_size[i], pRB_size[i],s.yStart,resultItemsIds);
                            }

                            // cout << "First " << resultItems.size() << endl;

                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], 0, pRA_size[i],s.xStart,resultItemsIds);
                                twoLevelPlus::window::RangeWithBinarySearch(pRXEnd[i], pRB_size[i], pRC_size[i], s.xStart,resultItemsIds);
                            }

                            // cout << "Second " << resultItems.size() << endl;

                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch2(pRXStart[i], 0, pRA_size[i],s.xEnd,resultItemsIds);
                            }

                            // cout << "Third " << resultItems.size() << endl;

                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevelPlus::window::RangeWithBinarySearch2(pRYStart[i], 0, pRA_size[i],s.yEnd,resultItemsIds);
                            }
                            // cout << "Fourth " << resultItems.size() << endl;
                    
                            break;
                    }
                    
                    timeBorders = tim.stop();
                    
                    tim.start();
                    unsigned long long newRes = 0;
                    for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                        newRes++;
                    }
                    timeCounting = tim.stop();


                    resultItemsIds.clear();
                    cout << "twoLevelPlus-Window-Query\t"<< j << "\t" << newRes <<"\t" << timeIndexingOrPartitioning <<"\t"<<timeInside+timeCorner+timeBorders+timeCounting +timeFindingTiles << endl;
                }
            }

            delete[] pRA_size;
            delete[] pRB_size;
            delete[] pRC_size;
            delete[] pRD_size;
            delete[] pRXStart;
            delete[] pRXEnd;
            delete[] pRYStart;
            delete[] pRYEnd;
            delete[] pR;
            
            break;

        case DISK_QUERY:
            cout << "Not supported!" << endl;
            break;

    }
    
}
