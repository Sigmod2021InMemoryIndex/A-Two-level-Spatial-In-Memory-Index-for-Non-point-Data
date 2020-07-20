#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./partitioning/partition.h"
#include "./grid/twoLevel.h"

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./twoLevel - range query using the 2-level algorithm" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./twoLevel [OPTION]... [FILE1] [FILE2]" << endl << endl;
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
    cerr << "       Window range query using the 2-level algorithm with 3000 partitions per dimension." << endl;
    cerr << "              ./twoLevel -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "       Disk range query using the 2-level algorithm with 3000 partitions and a radius of 0.1" << endl;
    cerr << "              ./twoLevel -p 3000 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    char c;
    int runNumPartitionsPerRelation = -1;
    Timer tim;
    double timeCounting = 0, timeIndexingOrPartitioning = 0 , timeFindingTiles = 0, timeCorner = 0,timeInside = 0,timeBorders = 0;
    Relation R, S, *pR;
    size_t *pRA_size, *pRB_size, *pRC_size, *pRD_size;
    int runNumPartitions = -1;
    int queryCase, queryMethod=-1;
    vector<RecordId> resultItemsIds;
    Coord epsilon = -1.0;
    Coord minX, maxX, minY, maxY, diffX, diffY, maxExtend;
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
                NUM_ITERATIONS = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cerr << "Wrong arguments! ";
                usage();
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

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;

    cout << "Method\t" << "Query ID\t" << "Results\t" << "Indexing Time\t" << "Filtering Time" << endl;

    switch (queryMethod)
    {
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

            pR = new Relation[runNumPartitions];

            pRA_size = new size_t[runNumPartitions];
            pRB_size = new size_t[runNumPartitions];
            pRC_size = new size_t[runNumPartitions];
            pRD_size = new size_t[runNumPartitions];
            memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRD_size, 0, runNumPartitions*sizeof(size_t));    
            
            tim.start();
            partition::twoLevel::single::PartitionTwoDimensional(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();

            for (size_t j = 0; j < S.numRecords; j++){
                for(int o = 0; o < NUM_ITERATIONS; o++) {
                    vector<int> insideCells;
                    vector<int> cornerCells;

                    result = 0;
                    auto &s = S[j];
                    
                    tim.start();
                    queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                    timeFindingTiles = tim.stop();
                    
                    auto lsize = insideCells.size();
                    for(int pid = 0; pid < lsize; pid++)
                    {
                        auto cid = insideCells[pid];
                        for (auto it = 0; it < pRA_size[cid]; it++) {
                            resultItemsIds.push_back(pR[cid][it].id);
                        }
                    }

                    for(int i = 0; i < 4; i++)
                    {
                        if(cornerCells[i] != -1)
                            twoLevel::window::Range_Corners(pR[cornerCells[i]],s, 0, pRA_size[cornerCells[i]],resultItemsIds);
                    }
                    for(int i = 0; i < 2; i++)
                    {
                        if(cornerCells[i] != -1)
                            twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRA_size[cornerCells[i]], pRB_size[cornerCells[i]],resultItemsIds);
                    }
                    for(int i = 0; i < 3; i+=2)
                    {
                        if(cornerCells[i] != -1)
                            twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRB_size[cornerCells[i]], pRC_size[cornerCells[i]],resultItemsIds);
                    }
                    twoLevel::window::Range_Corners(pR[cornerCells[0]],s, pRC_size[cornerCells[0]], pRD_size[cornerCells[0]],resultItemsIds);
                    
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds);
                            }
                            
                            break;
                        case 3:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds);
                            }
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds);
                            }
                            
                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                            }
                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                            }
                            
                            break;
                    }
                    
                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);

                    tim.start();
                    lsize = insideCells.size();
                    for(int pid = 0; pid < lsize; pid++)
                    {
                        auto cid = insideCells[pid];
                        for (auto it = 0; it < pRA_size[cid]; it++) {
                            resultItemsIds.push_back(pR[cid][it].id);
                        }
                    }
                    timeInside = tim.stop();

                    tim.start();
                    for(int i = 0; i < 4; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],resultItemsIds);
                    }
                    for(int i = 0; i < 2; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],resultItemsIds);
                    }
                    for(int i = 0; i < 3; i+=2)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1)
                            twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],resultItemsIds);
                    }
                    auto cid = cornerCells[0];
                    twoLevel::window::Range_Corners(pR[cid],s, pRC_size[cid], pRD_size[cid],resultItemsIds);
                    timeCorner = tim.stop();
                    
                    tim.start();
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds);
                            }
                            
                            break;
                        case 3:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds);
                            }
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                                twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds);
                            }
                            
                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                            }

                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],resultItemsIds);
                            }
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
                    cout << "twoLevel-Window-Query\t"<< j << "\t" << newRes <<"\t" << timeIndexingOrPartitioning << "\t" << timeInside+timeCorner+timeBorders+timeCounting +timeFindingTiles << endl;
                }
            }   
            
            delete[] pRA_size;
            delete[] pRB_size;
            delete[] pRC_size;
            delete[] pRD_size;
            delete[] pR;
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

            pR = new Relation[runNumPartitions];
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            pRA_size = new size_t[runNumPartitions];
            pRB_size = new size_t[runNumPartitions];
            pRC_size = new size_t[runNumPartitions];
            pRD_size = new size_t[runNumPartitions];
            memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRD_size, 0, runNumPartitions*sizeof(size_t));    
            
            tim.start();
            partition::twoLevel::single::PartitionTwoDimensional(R, pR,pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();

            for (size_t j = 0; j < S.numRecords; j++){
                for(int o = 0; o < NUM_ITERATIONS; o++) {
                    vector<int> insideCells;
                    vector<int> cornerCells;

                    result = 0;
                    auto &s = S[j];
                    
                    Coord xCircle,yCircle;

                    xCircle = (s.xStart + s.xEnd)/2;
                    yCircle = (s.yStart + s.yEnd)/2;
                    
                    Coordinates point (xCircle,yCircle);

                    tim.start();
                    queryCase = twoLevel::disk::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                    timeFindingTiles = tim.stop();

                    vector<RecordId> resultItemsIds;

                    auto lsize = insideCells.size();
                    for(int pid = 0; pid < lsize; pid++)
                    {
                        auto cid = insideCells[pid];

                        int x,y;
                        twoLevel::disk::findCoordinates(cid , x, y, runNumPartitionsPerRelation);

                        Coord xStart = x *partitionExtent;
                        Coord yStart = y *partitionExtent;

                        Record rec (0, xStart , yStart , xStart + partitionExtent, yStart + partitionExtent);
                        
                        if (twoLevel::disk::MaxDist(point, rec) <= epsilon*epsilon){
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                resultItemsIds.push_back(pR[cid][it].id);
                            }
                        }
                        else{
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                if (twoLevel::disk::MinDist(point,pR[cid][it]) > epsilon*epsilon){
                                    continue;
                                }
                                resultItemsIds.push_back(pR[cid][it].id);
                            }
                        }
                    }

                    for(int i = 0; i < 4; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){
                            twoLevel::disk::RangeCornersDisk(pR[cid],s, 0, pRA_size[cid],resultItemsIds,epsilon, point);
                        }
                    }

                    for(int i = 0; i < 2; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){

                            twoLevel::disk::RangeCornersDisk(pR[cid],s, pRA_size[cid], pRB_size[cid],resultItemsIds,epsilon, point);

                        }
                    }

                    for(int i = 0; i < 3; i+=2)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){

                            twoLevel::disk::RangeCornersDisk(pR[cid],s, pRB_size[cid], pRC_size[cid],resultItemsIds,epsilon, point);

                        }
                    }
                    auto cid = cornerCells[0];

                    twoLevel::disk::RangeCornersDisk(pR[cid],s, pRC_size[cid], pRD_size[cid],resultItemsIds,epsilon, point);

                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeBClassDisk(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            break;
                        case 2:
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeCClassDisk(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            break;
                        case 3:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeBClassDisk(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeCClassDisk(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                            }

                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                            }
                            break;
                    }

                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);

                    tim.start();
                    lsize = insideCells.size();
                    for(int pid = 0; pid < lsize; pid++)
                    {
                        auto cid = insideCells[pid];

                        int x,y;
                        twoLevel::disk::findCoordinates(cid , x, y, runNumPartitionsPerRelation);

                        Coord xStart = x *partitionExtent;
                        Coord yStart = y *partitionExtent;

                        Record rec (0, xStart , yStart , xStart + partitionExtent, yStart + partitionExtent);

                        if (twoLevel::disk::MaxDist(point, rec) <= epsilon*epsilon){
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                resultItemsIds.push_back(pR[cid][it].id);
                            }
                        }
                        else{
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                if (twoLevel::disk::MinDist(point,pR[cid][it]) > epsilon*epsilon){
                                    continue;
                                }
                                resultItemsIds.push_back(pR[cid][it].id);
                            }
                        }
                    }
                    timeInside = tim.stop();

                    
                    tim.start();
                    for(int i = 0; i < 4; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){
                            twoLevel::disk::RangeCornersDisk(pR[cid],s, 0, pRA_size[cid],resultItemsIds,epsilon, point);
                        }
                    }
                    for(int i = 0; i < 2; i++)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){
                            twoLevel::disk::RangeCornersDisk(pR[cid],s, pRA_size[cid], pRB_size[cid],resultItemsIds,epsilon, point);
                        }
                    }
                    for(int i = 0; i < 3; i+=2)
                    {
                        auto cid = cornerCells[i];
                        if(cid != -1){
                            twoLevel::disk::RangeCornersDisk(pR[cid],s, pRB_size[cid], pRC_size[cid],resultItemsIds,epsilon, point);
                        }
                    }
                    cid = cornerCells[0];

                    twoLevel::disk::RangeCornersDisk(pR[cid],s, pRC_size[cid], pRD_size[cid],resultItemsIds,epsilon, point);

                    timeCorner = tim.stop();
                    
                    tim.start();
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeBClassDisk(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds, epsilon, point);
                            }
                            break;
                        case 2:
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeCClassDisk(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            break;
                        case 3:
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeBClassDisk(pR[i],s, pRA_size[i], pRB_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                                twoLevel::disk::RangeCClassDisk(pR[i],s, pRB_size[i], pRC_size[i],resultItemsIds, epsilon, point);
                            }
                            
                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                twoLevel::disk::RangeCClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                            }

                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                twoLevel::disk::RangeBClassDisk(pR[i],s, 0, pRA_size[i],resultItemsIds, epsilon, point);
                            }
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
                    cout << "twoLevel-Disk-Query\t"<<j<<"\t" << "\t" << newRes <<"\t" << timeIndexingOrPartitioning <<"\t"<<timeInside+timeCorner+timeBorders+timeCounting +timeFindingTiles << endl;
                }
            }

            delete[] pRA_size;
            delete[] pRB_size;
            delete[] pRC_size;
            delete[] pRD_size;
            delete[] pR;
            
            break;
    } 
}
