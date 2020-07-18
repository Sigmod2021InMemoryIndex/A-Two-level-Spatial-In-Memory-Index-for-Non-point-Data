#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./grid/oneLevel.h"
#include "./partitioning/partition.h"

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./oneLevel - range query using the 1-level algorithm" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./oneLevel [OPTION]... [FILE1] [FILE2]" << endl << endl;
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
    cerr << "       Window range query using the 1-level algorithm with 3000 partitions per dimension." << endl;
    cerr << "              ./oneLevel -p 3000 -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "       Disk range query using the 1-level algorithm with 3000 partitions and a radius of 0.1" << endl;
    cerr << "              ./oneLevel -p 3000 -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    char c;
    int runNumPartitionsPerRelation = -1;
    Timer tim;
    double timeCounting = 0, timeIndexingOrPartitioning = 0, timeFindingTiles = 0, timeCorner = 0, timeInside = 0, timeBorders = 0,timeDisc = 0;
    Relation R, S, *pR;
    int runNumPartitions = -1;
    int queryCase, queryMethod = -1;
    vector<RecordId> resultItemsIds;
    Coord epsilon = -1.0;
    Coord minX, maxX, minY, maxY;
    Coord diffX, diffY, maxExtend;
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
            cerr<<"Radius value is missing"<<endl;
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
            
            tim.start();
            partition::oneLevel::single::PartitionTwoDimensional(R, pR, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();

            for (size_t j = 0; j < S.numRecords; j++){
                for(int o = 0; o < NUM_ITERATIONS; o++) {
                    auto &s = S[j];
                    result = 0;
                    vector<int> insideCells;
                    vector<int> cornerCells;


                    tim.start();
                    queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                    timeFindingTiles = tim.stop();

                    
                    for(int item = 0; item < insideCells.size(); item++)         
                    {
                        oneLevel::window::RangeInside(pR,insideCells[item],s,runNumPartitionsPerRelation,resultItemsIds);
                    }
                    
                    for(int item = 0; item < cornerCells.size() ; item++)         
                    {
                        if ( cornerCells[item] != -1){
                            oneLevel::window::RangeCorner(pR,cornerCells[item],s,runNumPartitionsPerRelation,resultItemsIds);
                        }
                    }

                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                oneLevel::window::RangeCorner(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeCorner(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                            
                            break;
                        case 3:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                oneLevel::window::RangeBottomBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }

                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeLeftBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                        

                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeRightBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                        
                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                oneLevel::window::RangeTopBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }

                            break;
                    }

                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);

                    tim.start();
                    auto lsize = insideCells.size();
                    for(int item = 0; item < lsize ; item++)
                    {
                        oneLevel::window::RangeInside(pR,insideCells[item],s,runNumPartitionsPerRelation,resultItemsIds);
                    }
                    timeInside = tim.stop();
                    
                    tim.start();
                    lsize = cornerCells.size();
                    for(int item = 0; item < lsize; item++)
                    {
                        auto cid = cornerCells[item];
                        if ( cid != -1){
                            oneLevel::window::RangeCorner(pR,cid,s,runNumPartitionsPerRelation,resultItemsIds);
                        }
                    }
                    timeCorner = tim.stop();

                    
                    tim.start();
                    switch (queryCase){
                        case 0: 
                            break;
                        case 1:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                oneLevel::window::RangeBottomBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                            
                            break;
                        case 2:
                            
                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeLeftBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                            
                            break;
                        case 3:
                            
                            for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                            {
                                oneLevel::window::RangeBottomBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }

                            for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeLeftBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                        

                            for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                            {
                                oneLevel::window::RangeRightBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
                            }
                        
                            for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                            {
                                oneLevel::window::RangeTopBorder(pR,i,s,runNumPartitionsPerRelation,resultItemsIds);
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
                    cout << "OneLevel-Window-Query\t"<< j << "\t" << newRes <<"\t" << timeIndexingOrPartitioning <<"\t"<<timeInside+timeCorner+timeBorders+timeCounting+timeFindingTiles <<  endl;
                }
            }
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

            pR = new Relation[runNumPartitions];
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            S.loadDisk(epsilon);

            tim.start();
            partition::oneLevel::single::PartitionTwoDimensional(R, pR, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();

            result = 0;
            
            for ( int i = 0 ; i < S.size(); i++){
                for(int o = 0; o < NUM_ITERATIONS; o++) {
                    auto &s = S[i];

                    Relation cellsCoord;
                    
                    tim.start();							
                    oneLevel::disk::FindRelevantTiles(s,cellsCoord, runNumPartitionsPerRelation);
                    timeFindingTiles = tim.stop();

                    
                    Coord xCircle,yCircle;
        
                    xCircle = (S[i].xStart + S[i].xEnd)/2;
                    yCircle = (S[i].yStart + S[i].yEnd)/2;
                    
                    Coordinates point2 (xCircle,yCircle);
                    
                    for ( int j = 0 ; j < cellsCoord.size() ; j ++ ){
                        for ( int k = 0 ; k < pR[cellsCoord[j].id].size(); k++){
                            auto &r = pR[cellsCoord[j].id][k];
                            if ( oneLevel::disk::MaxDist(point2,cellsCoord[j],2) <= epsilon*epsilon ){
                                auto x = max(s.xStart, r.xStart);
                                auto y = min(s.yEnd, r.yEnd);
                                
                                auto pid_ref = oneLevel::findReferenceCell(x, y, partitionExtent, runNumPartitionsPerRelation);
                                if (pid_ref == cellsCoord[j].id){
                                    resultItemsIds.push_back(s.id);
                                }
                            }
                            else{
                                if ((r.xStart > s.xEnd) || (r.xEnd < s.xStart) || (r.yStart > s.yEnd) || (r.yEnd < s.yStart)){
                                    continue;
                                }
                                else{
                                    auto x = max(s.xStart, r.xStart);
                                    auto y = min(s.yEnd, r.yEnd);
                                    
                                    auto pid_ref = oneLevel::findReferenceCell(x, y, partitionExtent, runNumPartitionsPerRelation);
                                    if (pid_ref == cellsCoord[j].id){

                                        if ( oneLevel::disk::MinDist(point2,r) <= epsilon*epsilon ){
                                            resultItemsIds.push_back(s.id);
                                        }
                                    }
                                }
                            }

                        }
                        
                    }

                    result = resultItemsIds.size();
                    resultItemsIds.clear(); 
                    resultItemsIds.reserve(result);
                    
                    tim.start();
                    xCircle = (S[i].xStart + S[i].xEnd)/2;
                    yCircle = (S[i].yStart + S[i].yEnd)/2;
                    
                    Coordinates point (xCircle,yCircle);
                    
                    for ( int j = 0 ; j < cellsCoord.size() ; j ++ ){
                        for ( int k = 0 ; k < pR[cellsCoord[j].id].size(); k++){
                            auto &r = pR[cellsCoord[j].id][k];
                            if ( oneLevel::disk::MaxDist(point,cellsCoord[j],2) <= epsilon*epsilon ){
                                auto x = max(s.xStart, r.xStart);
                                auto y = min(s.yEnd, r.yEnd);
                                
                                auto pid_ref = oneLevel::findReferenceCell(x, y, partitionExtent, runNumPartitionsPerRelation);
                                if (pid_ref == cellsCoord[j].id){
                                    resultItemsIds.push_back(s.id);
                                }
                            }
                            else{
                                if ((r.xStart > s.xEnd) || (r.xEnd < s.xStart) || (r.yStart > s.yEnd) || (r.yEnd < s.yStart)){
                                    continue;
                                }
                                else{
                                    auto x = max(s.xStart, r.xStart);
                                    auto y = min(s.yEnd, r.yEnd);
                                    
                                    auto pid_ref = oneLevel::findReferenceCell(x, y, partitionExtent, runNumPartitionsPerRelation);
                                    if (pid_ref == cellsCoord[j].id){

                                        if ( oneLevel::disk::MinDist(point,r) <= epsilon*epsilon ){
                                            resultItemsIds.push_back(s.id);
                                        }
                                    }
                                }
                            }

                        }
                        
                    }
                
                
                    timeDisc = tim.stop();
                        
                    tim.start();
                    unsigned long long newRes = 0;
                    for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                        newRes++;
                    }
                    timeCounting = tim.stop();

                    
                    resultItemsIds.clear();
                    cout << "oneLevel-Disk-Query\t"<< i << "\t" << newRes <<"\t" << timeIndexingOrPartitioning <<"\t"<<timeDisc+timeCounting+timeFindingTiles << endl;
                }
            }
            break;
    }
    
}
