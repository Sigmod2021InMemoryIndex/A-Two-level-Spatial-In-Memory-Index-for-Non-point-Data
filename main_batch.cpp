#include <cstdlib>
#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./grid/batchProcessing.h"
#include "./partitioning/partition.h"
using namespace std;


void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./batch - code for the batch processing experiments" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./batch [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       Mandatory arguments" << endl << endl;
    cerr << "       -p" << endl;
    cerr << "              number of partions per dimension" << endl;
    cerr << "       -t" << endl;
    cerr << "              number of threads" << endl;
    cerr << "       -v" << endl;
    cerr << "              batch processing method. 1 Queries-based, 2 for Tiles-based" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "        Queries-based batch processing window query using 1 thread." << endl;
    cerr << "              ./batch -p 3000 -t 1 -v 1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "        Tiles-based batch processing window query using 4 thread." << endl;
    cerr << "              ./batch -p 3000 -t 4 -v 2 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}


int main(int argc, char** argv) {
    char c;
    int runNumPartitionsPerRelation = -1;
    int runNumPartitions = -1;
    int runNumThreads = -1;
    Timer tim;
    double queryExecutionTime = 0,timeIndexingOrPartitioning = 0, timeFindingTiles = 0, timeTransform = 0;
    Relation R, S, *pRA, *pRB, *pRC, *pRD;
    Relation *pSIns, *pSCorDR, *pSCorDL, *pSCorUR, *pSCorUL, *pSBorX, *pSBorY, *pSBorB, *pSBorC ;
    size_t *pRA_size, *pRB_size, *pRC_size, *pRD_size;
    size_t *pSIns_size, *pSCorDR_size, *pSCorDL_size, *pSCorUR_size, *pSCorUL_size, *pSBorX_size, *pSBorY_size, *pSBorB_size, *pSBorC_size ;
    bool *checkTile;
    auto vsize = 0;
    vector<int> vecPart;
    int queryCase,parallelCase=-1;
    int NUM_ITERATIONS = 1;


    unsigned long long result = 0;

    while ((c = getopt(argc, argv, "p:t:v:i:m")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
                break;
            case 't':
                runNumThreads = atoi(optarg);
                break;
            case 'v':
                parallelCase = atoi(optarg);
                break;
            case 'i':
                NUM_ITERATIONS = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                break;
        }
    }

    if(runNumPartitionsPerRelation == -1)
    {
        cerr << "Number of partitions is missing" << endl;
        usage();
    }

    if(parallelCase == -1)
    {
        cerr << "Batch processing method is missing" << endl;
        usage();
    }

    if(runNumThreads == -1)
    {
        cerr << "Number of threads is missing" << endl;
        usage();
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
    
    
    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
    cout << "Method\t" << "# Threads\t" << "Total Results\t" << "Indexing Time\t" << "Total Filtering Time" << endl; 

    
    size_t results[S.numRecords];
    memset(results, 0, S.numRecords*sizeof(size_t));

    int queryCases[S.numRecords];
    vector<int> cells[S.numRecords];
    vector<int> insideCells[S.numRecords];
    vector<int> cornerCells[S.numRecords];


    for(int o = 0; o < NUM_ITERATIONS; o++) 
    {
        result = 0;
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

        if (runNumThreads == 1)
        {
            switch(parallelCase)
            {
                case QUERIES_BASED:
                    tim.start();
                    partition::update::PartitionTwoDimensional(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                    timeIndexingOrPartitioning = tim.stop();
                    
                    tim.start();
                    for (size_t j = 0; j < S.numRecords; j++){
                        vector<int> insideCells;
                        vector<int> cornerCells;
                        
                        result = 0;
                        auto &s = S[j];
                        
                        queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                        
                        auto lsize = insideCells.size();
                        for(int pid = 0; pid< lsize; pid++)
                        {
                            auto cid = insideCells[pid];
                            if ( pRA_size[insideCells[pid]] > 0 ){
                                for (auto it = 0; it < pRA_size[cid]; it++) {
                                    result += pRA[cid][it].id ^ s.id;
                                }
                            }
                        }
                        
                        for(int i = 0; i < 4; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRA_size[cornerCells[i]] > 0){
                                    result += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], s, pRA_size[cid]);
                                }
                            }
                        }
                        for(int i = 0; i < 2; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRB_size[cornerCells[i]] > 0){
                                    result += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], s, pRB_size[cid]);
                                }
                            }
                        }
                        for(int i = 0; i < 3; i+=2)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRC_size[cornerCells[i]] > 0){
                                    result += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], s, pRC_size[cid]);
                                }
                            }
                        }
                        
                        auto cid = cornerCells[0];
                        if (pRD_size[cid] > 0){
                            result += batchProcessing::window::xor_workload::RangeCorners(pRD[cid], s, pRD_size[cid]);
                        }
                        
                        switch (queryCase){
                            case 0:
                                break;
                            case 1:
                                lsize = cornerCells[1];
                                for(int i = cornerCells[0]+1; i < lsize; i++)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }
                                    
                                    if ( pRB_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeBClass(pRB[i], s, pRB_size[i]);
                                    }
                                }
                                
                                break;
                            case 2:
                                
                                lsize = cornerCells[2];
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }

                                    if ( pRC_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeCClass(pRC[i], s, pRC_size[i]);
                                    }
                                }
                                
                                break;
                            case 3:
                                
                                lsize = cornerCells[1];
                                for(int i = cornerCells[0]+1; i < lsize; i++)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }
                                    
                                    if ( pRB_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeBClass(pRB[i], s, pRB_size[i]);
                                    }
                                }
                                
                                lsize = cornerCells[2];
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0){
                                        result += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }
                                    
                                    if ( pRC_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeCClass(pRC[i], s, pRC_size[i]);
                                    }
                                }
                                
                                lsize = cornerCells[3];
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        result += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }
                                }

                                for ( int i = cornerCells[2] + 1; i< lsize; i++)
                                {
                                    if ( pRA_size[i] > 0){
                                        result += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }
                                }
                        }
                        
                        results[j] = result;
                    }
                    queryExecutionTime = tim.stop();
                    
                    result = 0;
                    for ( int  i = 0; i < S.numRecords; i++){
                        result += results[i];
                        // results[i] = 0;
                    }
                    memset(results, 0, S.numRecords*sizeof(size_t));
                    cout << "Bacth-Processing-Queries-Based\t" << "\t" << runNumThreads << "\t" << result << "\t" << timeIndexingOrPartitioning << "\t" << queryExecutionTime << endl;
                    break;
                    
                case TILES_BASED:
                    pSIns = new Relation[runNumPartitions];
                    pSCorDR = new Relation[runNumPartitions];
                    pSCorDL = new Relation[runNumPartitions];
                    pSCorUR = new Relation[runNumPartitions];
                    pSCorUL = new Relation[runNumPartitions];
                    pSBorX = new Relation[runNumPartitions];
                    pSBorY = new Relation[runNumPartitions];
                    pSBorB = new Relation[runNumPartitions];
                    pSBorC = new Relation[runNumPartitions];
                    
                    pSIns_size = new size_t[runNumPartitions];
                    pSCorDR_size = new size_t[runNumPartitions];
                    pSCorDL_size = new size_t[runNumPartitions];
                    pSCorUR_size = new size_t[runNumPartitions];
                    pSCorUL_size = new size_t[runNumPartitions];
                    pSBorX_size = new size_t[runNumPartitions];
                    pSBorY_size = new size_t[runNumPartitions];
                    pSBorB_size = new size_t[runNumPartitions];
                    pSBorC_size = new size_t[runNumPartitions];
                    
                    checkTile = new bool[runNumPartitions];
                    
                    memset(pSIns_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorDR_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorDL_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorUR_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorUL_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorX_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorY_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorB_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorC_size, 0, runNumPartitions*sizeof(size_t));
                    memset(checkTile, false, runNumPartitions*sizeof(bool));
                    
                    tim.start();
                    partition::update::PartitionTwoDimensional(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                    timeIndexingOrPartitioning = tim.stop();
                    
                    tim.start();
                    partition::batch::FindRelevantTiles(S, pSIns, pSCorDR, pSCorDL, pSCorUR, pSCorUL, pSBorX, pSBorY, pSBorB, pSBorC, pSIns_size, pSCorDR_size, pSCorDL_size, pSCorUR_size, pSCorUL_size, pSBorX_size, pSBorY_size, pSBorB_size, pSBorC_size, checkTile, runNumPartitionsPerRelation);
                    //

                    //tim.start();
                    for ( int i = 0 ; i < runNumPartitions ; i ++ ){
                        if ( checkTile[i] == true ){
                            vecPart.push_back(i);
                        }
                    }
                    //timeTransform = tim.stop();

                    //tim.start();
                    vsize = vecPart.size();
                    for (auto i = 0; i < vsize; i++) { 
                        auto cid = vecPart[i];
                        auto lsize = pSIns[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSIns[cid][j].id;
                            auto tmp = results[qid];
                            for ( int k = 0 ; k < pRA_size[cid] ; k ++){

                                tmp += pRA[cid][k].id ^ qid;

                            }
                            results[qid] = tmp;
                        }

                        lsize = pSCorDL[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSCorDL[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorDL[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], pSCorDL[cid][j], pRB_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], pSCorDL[cid][j], pRC_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRD[cid], pSCorDL[cid][j], pRD_size[cid]);
                            results[qid] = tmp;
                        }

                        lsize = pSCorDR[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSCorDR[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorDR[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], pSCorDR[cid][j], pRB_size[cid]);
                            results[qid] = tmp;
                        }

                        lsize = pSCorUL[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSCorUL[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorUL[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], pSCorUL[cid][j], pRC_size[cid]);
                            results[qid] = tmp;
                        }

                        lsize = pSCorUR[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            
                            results[pSCorUR[cid][j].id] += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorUR[cid][j], pRA_size[cid]);
                            
                        }

                        lsize = pSBorX[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSBorX[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeBClass(pRA[cid], pSBorX[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeBClass(pRB[cid], pSBorX[cid][j], pRB_size[cid]);
                            results[qid] = tmp;
                        }

                        lsize = pSBorB[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            results[pSBorB[cid][j].id] += batchProcessing::window::xor_workload::RangeBClass(pRA[cid], pSBorB[cid][j], pRA_size[cid]);
                        }

                        lsize = pSBorY[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSBorY[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCClass(pRA[cid], pSBorY[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCClass(pRC[cid], pSBorY[cid][j], pRC_size[cid]);
                            results[qid] = tmp;
                        }

                        lsize = pSBorC[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            results[pSBorC[cid][j].id] += batchProcessing::window::xor_workload::RangeCClass(pRA[cid], pSBorC[cid][j], pRA_size[cid]);
                        }

                    }
                    queryExecutionTime = tim.stop();
                    
                    for ( int  i = 0; i < S.numRecords; i++){
                        result += results[i];
                    }
                    memset(results, 0, S.numRecords*sizeof(size_t));
                    vecPart.clear();
                    
                    cout << "Batch-Processing-Tiles-Based\t" << runNumThreads << "\t" << result << "\t" << timeIndexingOrPartitioning << "\t" << queryExecutionTime<< endl;

                    delete[] pSIns;
                    delete[] pSCorDR;
                    delete[] pSCorDL;
                    delete[] pSCorUR;
                    delete[] pSCorUL;
                    delete[] pSBorX;
                    delete[] pSBorY;
                    delete[] pSBorB;
                    delete[] pSBorC;
                    delete[] pSIns_size;
                    delete[] pSCorDR_size;
                    delete[] pSCorDL_size;
                    delete[] pSCorUR_size;
                    delete[] pSCorUL_size;
                    delete[] pSBorX_size;
                    delete[] pSBorY_size;
                    delete[] pSBorB_size;
                    delete[] pSBorC_size;
                    delete[] checkTile;
                    
                    break;
            }

        }
        else
        {
            switch(parallelCase)
            {
                case QUERIES_BASED:
                    tim.start();
                    partition::update::PartitionTwoDimensional(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                    timeIndexingOrPartitioning = tim.stop();

                    tim.start();
                    #pragma omp parallel for schedule(dynamic)
                    for (size_t j = 0; j < S.numRecords; j++){
                        vector<int> insideCells;
                        vector<int> cornerCells;

                        result = 0;
                        auto &s = S[j];
                        
                        int queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);

                        auto lsize = insideCells.size();
                        for(int pid = 0; pid< lsize; pid++)
                        {
                            auto cid = insideCells[pid];
                            if ( pRA_size[cid] > 0 ){
                                for (auto it = 0; it < pRA_size[cid]; it++) {
                                    results[j] += pRA[cid][it].id ^ s.id;
                                }
                            }
                        }

                        for(int i = 0; i < 4; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRA_size[cid] > 0){
                                    results[j] += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], s, pRA_size[cid]);
                                }
                            }
                        }

                        for(int i = 0; i < 2; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRB_size[cid] > 0){
                                    results[j] += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], s, pRB_size[cid]);
                                }
                            }

                        }

                        for(int i = 0; i < 3; i+=2)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1){
                                if ( pRC_size[cid] > 0){
                                    results[j] += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], s, pRC_size[cid]);
                                }
                            }
                        }

                        auto cid = cornerCells[0];
                        if (pRD_size[cid] > 0){
                            results[j] += batchProcessing::window::xor_workload::RangeCorners(pRD[cid], s, pRD_size[cid]);
                        }

                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                lsize = cornerCells[1];
                                for(int i = cornerCells[0]+1; i < lsize; i++)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }

                                    if ( pRB_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeBClass(pRB[i], s, pRB_size[i]);
                                    }
                                }

                                break;
                            case 2:
                                lsize = cornerCells[2];
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }
                                    if ( pRC_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeCClass(pRC[i], s, pRC_size[i]);
                                    }
                                }

                                break;
                            case 3:
                                lsize = cornerCells[1];
                                for(int i = cornerCells[0]+1; i < lsize; i++)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }

                                    if ( pRB_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeBClass(pRB[i], s, pRB_size[i]);
                                    }
                                }

                                lsize = cornerCells[2];
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0){
                                        results[j] += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }

                                    if ( pRC_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeCClass(pRC[i], s, pRC_size[i]);
                                    }
                                }

                                lsize = cornerCells[3];
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < lsize; i+=runNumPartitionsPerRelation)
                                {
                                    if ( pRA_size[i] > 0 ){
                                        results[j] += batchProcessing::window::xor_workload::RangeCClass(pRA[i], s, pRA_size[i]);
                                    }
                                }
                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    if ( pRA_size[i] > 0){
                                        results[j] += batchProcessing::window::xor_workload::RangeBClass(pRA[i], s, pRA_size[i]);
                                    }

                                }
                                break;
                        }
                    }
                    queryExecutionTime = tim.stop();
                    for ( int  i = 0; i < S.numRecords; i++){
                        result += results[i];
                    }
                    memset(results, 0, S.numRecords*sizeof(size_t));
                    cout << "Bacth-Processing-Queries-Based\t" << runNumThreads << "\t" << result <<"\t" << timeIndexingOrPartitioning <<"\t"<<queryExecutionTime<< endl;
                    break;
                    
                case TILES_BASED:
                    tim.start();
                    partition::update::PartitionTwoDimensional(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                    timeIndexingOrPartitioning = tim.stop();
                    
                    pSIns = new Relation[runNumPartitions];
                    pSCorDR = new Relation[runNumPartitions];
                    pSCorDL = new Relation[runNumPartitions];
                    pSCorUR = new Relation[runNumPartitions];
                    pSCorUL = new Relation[runNumPartitions];
                    pSBorX = new Relation[runNumPartitions];
                    pSBorY = new Relation[runNumPartitions];
                    pSBorB = new Relation[runNumPartitions];
                    pSBorC = new Relation[runNumPartitions];
                    
                    pSIns_size = new size_t[runNumPartitions];
                    pSCorDR_size = new size_t[runNumPartitions];
                    pSCorDL_size = new size_t[runNumPartitions];
                    pSCorUR_size = new size_t[runNumPartitions];
                    pSCorUL_size = new size_t[runNumPartitions];
                    pSBorX_size = new size_t[runNumPartitions];
                    pSBorY_size = new size_t[runNumPartitions];
                    pSBorB_size = new size_t[runNumPartitions];
                    pSBorC_size = new size_t[runNumPartitions];
                
                    checkTile = new bool[runNumPartitions];

                    memset(pSIns_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorDR_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorDL_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSCorUR_size, 0, runNumPartitions*sizeof(size_t)); 
                    memset(pSCorUL_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorX_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorY_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorB_size, 0, runNumPartitions*sizeof(size_t));
                    memset(pSBorC_size, 0, runNumPartitions*sizeof(size_t));
                    memset(checkTile, false, runNumPartitions*sizeof(bool));         

                    tim.start();
                    partition::batch::FindRelevantTiles(S, pSIns, pSCorDR, pSCorDL, pSCorUR, pSCorUL, pSBorX, pSBorY, pSBorB, pSBorC, pSIns_size, pSCorDR_size, pSCorDL_size, pSCorUR_size, pSCorUL_size, pSBorX_size, pSBorY_size, pSBorB_size, pSBorC_size, checkTile, runNumPartitionsPerRelation);
                    //timeFindingTiles = tim.stop();
                    
                    //tim.start();
                    for ( int i = 0 ; i < runNumPartitions ; i ++ ){
                        if ( checkTile[i] == true ){
                            vecPart.push_back(i);
                        }
                    }
                    //timeTransform = tim.stop();
                    
                    //tim.start();
                    vsize = vecPart.size();
                    #pragma omp parallel for reduction(+: results[:S.numRecords]) schedule(dynamic)
                    for (auto i = 0; i < vsize; i++) {
                        auto cid = vecPart[i];
                        auto lsize = pSIns[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSIns[cid][j].id;
                            auto tmp = results[qid];
                            for ( int k = 0 ; k < pRA_size[cid] ; k ++){
                                
                                tmp += pRA[cid][k].id ^ qid;
                                
                            }
                            results[qid] = tmp;
                        }
                        
                        lsize = pSCorDL[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSCorDL[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorDL[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], pSCorDL[cid][j], pRB_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], pSCorDL[cid][j], pRC_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRD[cid], pSCorDL[cid][j], pRD_size[cid]);
                            results[qid] = tmp;
                        }
                        
                        lsize = pSCorDR[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            
                            auto qid = pSCorDR[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorDR[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRB[cid], pSCorDR[cid][j], pRB_size[cid]);
                            results[qid] = tmp;
                        }
                        
                        lsize = pSCorUL[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSCorUL[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorUL[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCorners(pRC[cid], pSCorUL[cid][j], pRC_size[cid]);
                            results[qid] = tmp;
                        }
                        
                        lsize = pSCorUR[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            
                            results[pSCorUR[cid][j].id] += batchProcessing::window::xor_workload::RangeCorners(pRA[cid], pSCorUR[cid][j], pRA_size[cid]);
                            
                        }

                        lsize = pSBorX[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSBorX[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeBClass(pRA[cid], pSBorX[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeBClass(pRB[cid], pSBorX[cid][j], pRB_size[cid]);
                            results[qid] = tmp;
                        }
                        
                        lsize = pSBorB[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            results[pSBorB[cid][j].id] += batchProcessing::window::xor_workload::RangeBClass(pRA[cid], pSBorB[cid][j], pRA_size[cid]);
                        }
                        
                        lsize = pSBorY[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            auto qid = pSBorY[cid][j].id;
                            auto tmp = results[qid];
                            tmp += batchProcessing::window::xor_workload::RangeCClass(pRA[cid], pSBorY[cid][j], pRA_size[cid]);
                            tmp += batchProcessing::window::xor_workload::RangeCClass(pRC[cid], pSBorY[cid][j], pRC_size[cid]);
                            results[qid] = tmp;
                        }
                        
                        lsize = pSBorC[cid].size();
                        for ( int j = 0 ; j < lsize; j++ ){
                            results[pSBorC[cid][j].id] += batchProcessing::window::xor_workload::RangeCClass(pRA[cid], pSBorC[cid][j], pRA_size[cid]);
                        }
                    }
                    queryExecutionTime = tim.stop();
                    
                    
                    for ( int  i = 0; i < S.numRecords; i++){
                        result += results[i];
                    }
                    memset(results, 0, S.numRecords*sizeof(size_t));
                    vecPart.clear();

                    cout << "Bacth-Processing-Tiles-Based\t" << runNumThreads << "\t" << result << "\t" << timeIndexingOrPartitioning << "\t" << queryExecutionTime<< endl;
                    
                    delete[] pSIns;
                    delete[] pSCorDR;
                    delete[] pSCorDL;
                    delete[] pSCorUR;
                    delete[] pSCorUL;
                    delete[] pSBorX;
                    delete[] pSBorY;
                    delete[] pSBorB;
                    delete[] pSBorC;
                    delete[] pSIns_size;
                    delete[] pSCorDR_size;
                    delete[] pSCorDL_size;
                    delete[] pSCorUR_size;
                    delete[] pSCorUL_size;
                    delete[] pSBorX_size;
                    delete[] pSBorY_size;
                    delete[] pSBorB_size;
                    delete[] pSBorC_size;
                    delete[] checkTile;
                    
                    break;
            }
        }
        
        delete[] pRA;
        delete[] pRB;
        delete[] pRC;
        delete[] pRD;
        delete[] pRA_size;
        delete[] pRB_size;
        delete[] pRC_size;
        delete[] pRD_size;
    }
    return 0;
}


