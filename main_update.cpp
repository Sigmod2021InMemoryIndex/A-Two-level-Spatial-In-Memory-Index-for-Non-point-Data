#include <cstdlib>
#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./partitioning/partition.h"
#include "./grid/twoLevel.h"
#include "./grid/oneLevel.h"
#include "./quadTree/QuadTree.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>

#define MAX_NODE_CAPACITY   16

using namespace std;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./update - code for the update experiments" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./update [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       Mandatory arguments" << endl << endl;
    cerr << "       -q" << endl;
    cerr << "              updates using Quad-Tree. Option -c should be used when using -q" << endl;
    cerr << "       -c" << endl;
    cerr << "              node capacity" << endl;
    cerr << "       -1" << endl;
    cerr << "              updates using 1-Level. Option -p should be used when using -1 " << endl;
    cerr << "       -2" << endl;
    cerr << "              updates using 2-Level. Option -p should be used when using -2" << endl;
    cerr << "       -p" << endl;
    cerr << "              number of partitions per dimension" << endl;
    cerr << "       -r" << endl;
    cerr << "              updates using R-Tree " << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -h" << endl;
    cerr << "              heigh of Quad-Tree" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "        Updates using Quad-Tree with capacity 1000." << endl;
    cerr << "              ./update -c 1000 -q dataFile_90%.inp TIGER_c0.1%_n10000.qry dataFile_10%.inp" << endl;
    cerr << "        Updates using 2-Level with 3000 partitions per dimension" << endl;
    cerr << "              ./update -p 3000 -2 dataFile_90%.inp TIGER_c0.1%_n10000.qry dataFile_10%.inp" << endl;
    cerr << "\n" << endl;
    exit(1);
}




int main(int argc, char** argv) {
    char c;
    Relation R, S, U, *pRA, *pRB, *pRC, *pRD, *pR;
    size_t *pRA_size, *pRB_size, *pRC_size, *pRD_size, *pR_size;
    double timeQuery = 0;
    Timer tim;
    int algoMethod = -1;
    int runNumPartitions = -1;
    int runNumPartitionsPerRelation = -1;
    double timeIndexingOrPartitioning = 0, timeJoining = 0 , timeFindingTiles = 0, queryExecutionTime = 0, timeCounting = 0, timeInsert = 0;
    double timeCorner = 0,timeInside = 0,timeBorders = 0;
    int queryCase;
    int sizeOfR = 0;
    int sizeOfRAtFirst = 0;
    int NUM_ITERATIONS = 1; 

    typedef bg::model::point<double, 2, bg::cs::cartesian> point;
    typedef bg::model::box<point> box;
    typedef std::pair<box, unsigned> value;
    
    double totalTime;
    int numOfBoolVector = 0;
    unsigned long long result = 0;
    unsigned long long results = 0;
    unsigned long long totalRes = 0;
    vector<RecordId> resultItemsIds;
    int level = -1;
    QuadTreeNode* tree;

    

    int capacity = -1;

    while ((c = getopt(argc, argv, "p:rqc:h:12m")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
                break;
            case 'r':
                algoMethod = R_TREE;
                break;
            case 'q':
                algoMethod = QUAD_TREE;
                break;
            case '1':
                algoMethod = ONE_LEVEL;
                break;
            case '2':
                algoMethod = SECOND_LEVEL;
                break;
            case 'c':
                capacity = atoi(optarg);
                break;
            case 'h':
                level = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cout << "Wrong arguments! ";
                usage();
                break;
        }
    }

    if(algoMethod == -1)
    {
        cerr << "Query method is missing" << endl;
        usage();
    }

    if(algoMethod == ONE_LEVEL || algoMethod == SECOND_LEVEL)
    {
        if(runNumPartitionsPerRelation == -1)
        {
            cerr << "Number of partitions is missing" << endl;
            usage();
        }
    }

    if (algoMethod == QUAD_TREE){
        if(capacity == -1)
        {
            cout<<"Capacity value is missing"<<endl;
            usage();
        }
        if (capacity < 100){
            cout<<"Capacity value should be greater than 100"<<endl;
            usage();
        }
    }

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
    
    cout << "Method\t"<< "Results\t" << "Indexing Time\t" << "Update Query Time\t" <<"Window Query Time\t"<<"Total Time"<< endl;
    timeIndexingOrPartitioning = 0; 
    timeQuery = 0;
    timeInsert = 0;
    timeFindingTiles = 0;
        
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
        #pragma omp section
        {
            U.load(argv[optind+2]);
        }
    }

    Coord minX = min(R.minX, S.minX);
    Coord maxX = max(R.maxX, S.maxX);
    Coord minY = min(R.minY, S.minY);
    Coord maxY = max(R.maxY, S.maxY);
    
    minX = min(minX, U.minX);
    maxX = max(maxX, U.maxX);
    minY = min(minY, U.minY);
    maxY = max(maxY, U.maxY);
    
    Coord diffX = maxX - minX;
    Coord diffY = maxY - minY;
    Coord maxExtent = (diffX<diffY)?diffY:diffX;

    R.normalize(minX, maxX, minY, maxY, maxExtent);
    S.normalize(minX, maxX, minY, maxY, maxExtent);
    U.normalize(minX, maxX, minY, maxY, maxExtent);

    sizeOfRAtFirst = R.numRecords;
    switch (algoMethod){
        case R_TREE:{
            std::vector<box> contourCenters; 
            std::vector<value> cloud;

            for ( unsigned i = 0 ; i < R.size() ; ++i )
            {
                box b(point(R[i].xStart, R[i].yStart), point(R[i].xEnd, R[i].yEnd));
                contourCenters.push_back(b);
            }

            tim.start();
            size_t id_gen = 0;
            std::transform(
                    contourCenters.begin(), contourCenters.end(),
                    back_inserter(cloud), 
                    [&](box const& p) { return std::make_pair(p, id_gen++); }
            );

            bgi::rtree<value, bgi::rstar<MAX_NODE_CAPACITY> > rtree(cloud.begin(),cloud.end());

            timeIndexingOrPartitioning = tim.stop(); 

            vector<box> queries;
            for ( unsigned i = 0 ; i < S.size() ; ++i )
            {
                box query_box(point(S[i].xStart, S[i].yStart), point(S[i].xEnd, S[i].yEnd));
                queries.push_back(query_box);
            }

            int size = R.numRecords;
            
            tim.start();
            for (int i = 0 ; i < U.numRecords ; i ++ )
            {
                box b(point(U[i].xStart, U[i].yStart), point(U[i].xEnd, U[i].yEnd));
                rtree.insert(std::make_pair(b, size));
                size++;
            }
            timeInsert = tim.stop();
            
            vector<double> indexTime[NUM_ITERATIONS];
            vector<double> queryTime[NUM_ITERATIONS];
            for ( int i = 0 ; i < queries.size() ; i ++){
                unsigned long long results = 0;

                tim.start();

                std::vector<value> result_s;
                rtree.query(bgi::intersects(queries.at(i)), std::back_inserter(result_s));

                BOOST_FOREACH(value const& v, result_s)
                    results++;

                timeQuery += tim.stop();
                result += results;
            }
            
            cout << "r-tree-UPDATE\t"<< "\t" << result <<"\t" << timeIndexingOrPartitioning<<"\t"<<timeInsert<<"\t"<<timeQuery <<"\t"<<timeInsert+timeQuery << endl;
    
            break;
        }
        
        case QUAD_TREE:{
            BoundingBox boundary(0.0, 0.0,1.0,1.0);
            if( level == -1 ){
                tim.start();
                tree = new QuadTreeNode(boundary,capacity, 0);
                tree->loadDataNoHeight(R);
                timeIndexingOrPartitioning = tim.stop();

                sizeOfR = R.numRecords;

                //put dataset percentage in the R
                for ( int i = 0 ; i < U.numRecords; i ++){
                    R.emplace_back(sizeOfR, U[i].xStart , U[i].yStart , U[i].xEnd, U[i].yEnd);
                    R.numRecords++;
                    sizeOfR ++;
                }

                tim.start();
                for (int i = sizeOfRAtFirst ; i < sizeOfR ; i ++ )
                {   
                    auto &r = R[i];
                    tree->insertNoHeight(r);
                }
                timeInsert = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    auto &s = S[i];

                    tree->rangeQuery(s,resultItemsIds);

                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);
                    
                    tim.start();
                    tree->rangeQuery(s,resultItemsIds);
                   
                    unsigned long long newRes = 0;
                    
                    for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                        newRes++;
                    }
                    timeQuery += tim.stop();
                    
                    results += newRes;
                    resultItemsIds.clear();    
                }

            } 
            else
            {
                tim.start();
                tree = new QuadTreeNode(boundary, capacity, 0, level);
                tree->loadData(R);
                timeIndexingOrPartitioning = tim.stop();

                sizeOfR = R.numRecords;

                //put dataset percentage in the R
                for ( int i = 0 ; i < U.numRecords; i ++){
                    R.emplace_back(sizeOfR, U[i].xStart , U[i].yStart , U[i].xEnd, U[i].yEnd);
                    R.numRecords++;
                    sizeOfR ++;
                }

                tim.start();
                for (int i = sizeOfRAtFirst ; i < sizeOfR ; i ++ )
                {   
                    auto &r = R[i];
                    tree->insert(r);
                }
                timeInsert = tim.stop();

                for (size_t i = 0; i < S.numRecords; i++){
                    auto &s = S[i];

                    tree->rangeQuery(s, resultItemsIds);

                    result = resultItemsIds.size();
                    resultItemsIds.clear();
                    resultItemsIds.reserve(result);
                    
                    tim.start();
                    tree->rangeQuery(s,resultItemsIds);
                    
                    unsigned long long newRes = 0;

                    for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                        newRes++;
                    }
                    timeQuery += tim.stop();

                    results += newRes;
                    resultItemsIds.clear();
                    
                }
            }
            if( level == -1 ){
                cout << "QuadTreeNoHeight-UPDATE\t" << results <<"\t" << timeIndexingOrPartitioning<<"\t"<<timeInsert<<"\t"<<timeQuery <<"\t"<<timeQuery+timeInsert << endl;
            }
            else{
                cout << "QuadTree-UPDATE\t" << results <<"\t" << timeIndexingOrPartitioning<<"\t"<<timeInsert<<"\t"<<timeQuery <<"\t"<<timeQuery+timeInsert << endl;
            }

            break;
        }

        case ONE_LEVEL:
            pR = new Relation[runNumPartitions];

            pR_size = new size_t[runNumPartitions];
            memset(pR_size, 0, runNumPartitions*sizeof(size_t));

            //load dataset percentage and calculate resize of partitions
            for (int i = 0 ; i < U.numRecords ; i ++ )
            {
               partition::update::UpdateCounters(U[i], pR, pR_size, runNumPartitionsPerRelation);
            }

            //partition of R
            tim.start();
            partition::update::PartitionTwoDimensional(R, pR, pR_size, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();
            
            sizeOfR = R.numRecords;

            //put dataset percentage in the R
            for ( int i = 0 ; i < U.numRecords; i ++){
                R.emplace_back(sizeOfR, U[i].xStart , U[i].yStart , U[i].xEnd, U[i].yEnd);
                R.numRecords++;
                sizeOfR ++;
            }

            //insert dataset percentage in the partitions
            tim.start();
            for (int i = sizeOfRAtFirst ; i < sizeOfR ; i ++ )
            {
                partition::update::insert(R[i], pR, pR_size, runNumPartitionsPerRelation);
            }
            timeInsert = tim.stop();

            //execute queries
            for (size_t j = 0; j < S.numRecords; j++){
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

                unsigned long long newRes = 0;
                for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                    newRes++;
                }
                timeQuery += tim.stop();

                results += newRes;
                resultItemsIds.clear();
            }
            cout << "oneLevel-Update\t"<< "\t" << results <<"\t" << timeIndexingOrPartitioning<<"\t"<<timeInsert<<"\t"<< timeQuery<<"\t"<<timeQuery+timeInsert+timeFindingTiles << endl;

            break;
        case SECOND_LEVEL: 
        
            pRA = new Relation[runNumPartitions];
            pRB = new Relation[runNumPartitions];
            pRC = new Relation[runNumPartitions];
            pRD = new Relation[runNumPartitions];

            pRA_size = new size_t[runNumPartitions];
            pRB_size = new size_t[runNumPartitions];
            pRC_size = new size_t[runNumPartitions];
            pRD_size = new size_t[runNumPartitions];
            memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
            memset(pRD_size, 0, runNumPartitions*sizeof(size_t)); 

            //load dataset percentage and calculate resize of partitions
            for (int i = 0 ; i < U.numRecords ; i ++ )
            {
                partition::update::UpdateCounters(U[i], pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
            }

            //partition of R
            tim.start();
            partition::update::PartitionTwoDimensional(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
            timeIndexingOrPartitioning = tim.stop();
            
            sizeOfR = R.numRecords;

            //put dataset percentage in the R
            for ( int i = 0 ; i < U.numRecords; i ++){
                R.emplace_back(sizeOfR, U[i].xStart , U[i].yStart , U[i].xEnd, U[i].yEnd);
                sizeOfR ++;
                R.numRecords ++;
            }

            //insert dataset percentage in the partitions
            tim.start();
            for (int i = sizeOfRAtFirst ; i < sizeOfR ; i ++ )
            {
                partition::update::insert(R[i], pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
            }
            timeInsert = tim.stop();
            
            //execute queries
            for (size_t j = 0; j < S.numRecords; j++){
                vector<int> insideCells;
                vector<int> cornerCells;

                auto &s = S[j];

                tim.start();
                queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                timeFindingTiles += tim.stop();
                
                for(int pid = 0; pid < insideCells.size(); pid++)
                {
                    auto cid = insideCells[pid];
                    if ( pRA_size[cid] > 0 ){
                        for (auto it = 0; it < pRA_size[cid]; it++) {
                            resultItemsIds.push_back(pRA[cid][it].id);
                        }
                    }
                }

                for(int i = 0; i < 4; i++)
                {
                    if(cornerCells[i] != -1){
                        if( pRA_size[cornerCells[i]] > 0 ){
                            twoLevel::update::Range_Corners(pRA[i],s, pRA_size[i],resultItemsIds);
                        }
                    }
                }
                
                for(int i = 0; i < 2; i++)
                {
                    if(cornerCells[i] != -1){
                        if( pRB_size[cornerCells[i]] > 0 ){
                            twoLevel::update::Range_Corners(pRB[i],s, pRB_size[i],resultItemsIds);
                        }
                    }
                }

                for(int i = 0; i < 3; i+=2)
                {
                    if(cornerCells[i] != -1){
                        if( pRC_size[cornerCells[i]] > 0 ){
                            twoLevel::update::Range_Corners(pRC[i],s, pRC_size[i],resultItemsIds);
                        }
                    }
                }


                if( pRD_size[cornerCells[0]] > 0 ){
                    twoLevel::update::Range_Corners(pRD[cornerCells[0]],s, pRD_size[cornerCells[0]],resultItemsIds);
                }

                switch (queryCase){
                    case 0:
                        break;
                    case 1:
                        for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if ( pRB_size[i] > 0 ){
                                twoLevel::update::Range_B_Class(pRB[i],s,  pRB_size[i],resultItemsIds);
                            }
                        }
                        break;
                    case 2:
                        for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if ( pRC_size[i] > 0 ){
                                twoLevel::update::Range_C_Class(pRC[i],s,  pRC_size[i],resultItemsIds);
                            }
                        }
                        break;
                    case 3:
                        for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if ( pRB_size[i] > 0 ){
                                twoLevel::update::Range_B_Class(pRB[i],s, pRB_size[i],resultItemsIds);
                            }
                        }

                        for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if ( pRC_size[i] > 0 ){
                                twoLevel::update::Range_C_Class(pRC[i],s, pRC_size[i],resultItemsIds);
                            }
                        }

                        for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }
                        }

                        for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                        {
                            if ( pRA_size[i] > 0 ){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                        }
                        break;
                }


                result = resultItemsIds.size();
                resultItemsIds.clear();
                resultItemsIds.reserve(result);

                tim.start();
                auto lsize = insideCells.size();
                for(int pid = 0; pid < lsize; pid++)
                {
                    auto cid = insideCells[pid];
                    if (pRA_size[cid] > 0){
                        for (auto it = 0; it < pRA_size[cid]; it++) {
                            resultItemsIds.push_back(pRA[cid][it].id);
                        }
                    }
                }

                for(int i = 0; i < 4; i++)
                {
                    auto cid = cornerCells[i];
                    if(cid != -1){
                        if (pRA_size[cid] > 0){
                            twoLevel::update::Range_Corners(pRA[cid],s, pRA_size[cid],resultItemsIds);
                        }
                    }
                }

                for(int i = 0; i < 2; i++)
                {
                    auto cid = cornerCells[i];
                    if(cid != -1){
                        if (pRB_size[cid] > 0){
                            twoLevel::update::Range_Corners(pRB[cid],s, pRB_size[cid],resultItemsIds);
                        }
                    }
                }

                for(int i = 0; i < 3; i+=2)
                {
                    auto cid = cornerCells[i];
                    if(cid != -1){
                        if (pRC_size[cid] > 0){
                            twoLevel::update::Range_Corners(pRC[cid],s, pRC_size[cid],resultItemsIds);
                        }
                    }
                }

                auto cid = cornerCells[0];
                if (pRD_size[cid] > 0){
                    twoLevel::update::Range_Corners(pRD[cid],s, pRD_size[cid],resultItemsIds);
                }

                switch (queryCase){
                    case 0:
                        break;
                    case 1:
                        for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if (pRB_size[i] > 0){
                                twoLevel::update::Range_B_Class(pRB[i],s,  pRB_size[i],resultItemsIds);
                            }
                        }
                        break;
                    case 2:
                        for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if (pRC_size[i] > 0){
                                twoLevel::update::Range_C_Class(pRC[i],s,  pRC_size[i],resultItemsIds);
                            }
                        }
                        break;
                    case 3:
                        for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if (pRB_size[i] > 0){
                                twoLevel::update::Range_B_Class(pRB[i],s, pRB_size[i],resultItemsIds);
                            }
                        }

                        for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }

                            if (pRC_size[i] > 0){
                                twoLevel::update::Range_C_Class(pRC[i],s, pRC_size[i],resultItemsIds);
                            }
                        }

                        for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_C_Class(pRA[i],s, pRA_size[i],resultItemsIds);
                            }
                        }

                        for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                        {
                            if (pRA_size[i] > 0){
                                twoLevel::update::Range_B_Class(pRA[i],s, pRA_size[i],resultItemsIds);                                }

                        }
                        break;
                }

                unsigned long long newRes = 0;
                for (auto it = resultItemsIds.begin(); it != resultItemsIds.end(); ++it) {
                    newRes++;
                }
                timeQuery += tim.stop();

                totalRes += newRes;
                resultItemsIds.clear();
            }
            cout << "twoLevel-Update\t"<< "\t" << totalRes <<"\t" << timeIndexingOrPartitioning<<"\t"<<timeInsert<<"\t"<<timeQuery<<"\t"<<timeQuery+timeInsert+timeFindingTiles<< endl;
        break; 
    }
        
    return 0;
}
