#include "def.h"
#include "getopt.h"
#include "./containers/relation.h"
#include "./partitioning/partition.h"
#include "./grid/twoLevel.h"
#include "./grid/refinement.h"

void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./refine - code for the refinement experiments" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./refine [OPTION]... [FILE1] [FILE2]" << endl << endl;
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
    cerr << "       -v" << endl;
    cerr << "              refinement method. 1 for Simple, 2 for RefAvoid, 3 for RefAvoid+" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "       Refinement on a window query using the RefAvoid method and 3000 partitions per dimension." << endl;
    cerr << "              ./refine -p 3000 -w -v 2 TIGER_ROADS_geom.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "        Refinement on a window query using the RefAvoid method and 3000 partitions and a radius of 0.1" << endl;
    cerr << "              ./refine -p 3000 -d -e 0.1 -v 2 TIGER_ROADS_geom.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}



int main(int argc, char **argv)
{
    char c;
    int runAlgorithm = -1;
    int runNumPartitionsPerRelation = -1;
    Timer tim;
    double timeRefinement = 0, timeCounting = 0, timeIndexingOrPartitioning = 0,timeFindingTiles = 0, timeCorner = 0,timeInside = 0,timeBorders = 0;
    Relation R, S, *pR;
    size_t *pRA_size, *pRB_size, *pRC_size, *pRD_size;
    int runNumPartitions = -1;
    int queryCase, refinementCase=-1, queryMethod=-1;
    Coord epsilon = -1.0;
    Coord minX, maxX, minY, maxY, diffX, diffY, maxExtend;
    int NUM_ITERATIONS = 1;


    unsigned long long result = 0;

    while ((c = getopt(argc, argv, "p:wdue:v:i:m")) != -1)
    {
        switch (c)
        {
            case 'p':
                runNumPartitionsPerRelation = atoi(optarg);
                break;
            case 'w':
                queryMethod = WINDOW_QUERY;
                break;
            case 'd':
                queryMethod = DISK_QUERY;
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'v':
                refinementCase = atoi(optarg);
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

    if(refinementCase == -1)
    {
        cerr << "Refinement method is missing" << endl;
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

    vector<Geometry> geo;
    #pragma omp parallel sections
    {
         #pragma omp section
         {
            R.loadWithGeometries(argv[optind],geo);
         }
         #pragma omp section
         {
            S.load(argv[optind+1]);
         }
    }

    cout << "Method\t" << "Query ID\t" << "Results\t" << "Results After Refinement\t" << "Items Need Refinement\t" << "Filtering Time\t" << "Secondary Filtering Time\t" << "Refinement Time\t"<<"Total Time" << endl;

    runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
    bool flag = false;
    
    int size = S.numRecords;

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
        refinement::normalizeGeometries(minX,minY,maxExtend,geo);

        pR = new Relation[runNumPartitions];

        pRA_size = new size_t[runNumPartitions];
        pRB_size = new size_t[runNumPartitions];
        pRC_size = new size_t[runNumPartitions];
        pRD_size = new size_t[runNumPartitions];
        memset(pRA_size, 0, runNumPartitions*sizeof(size_t));
        memset(pRB_size, 0, runNumPartitions*sizeof(size_t));
        memset(pRC_size, 0, runNumPartitions*sizeof(size_t));
        memset(pRD_size, 0, runNumPartitions*sizeof(size_t));    
        
        switch(refinementCase)
        {
            case SIMPLE_REFINEMENT:
                tim.start();
                partition::twoLevel::single::PartitionTwoDimensional(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t j = 0; j < size; j++)
                {
                    for(int i = 0; i < NUM_ITERATIONS; i++)
                    {
                        vector<RecordId> candItems;
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
                                candItems.push_back(pR[cid][it].id);
                            }
                        }

                        for(int i = 0; i < 4; i++)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, 0, pRA_size[cornerCells[i]],candItems);
                        }
                        for(int i = 0; i < 2; i++)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRA_size[cornerCells[i]], pRB_size[cornerCells[i]],candItems);
                        }
                        for(int i = 0; i < 3; i+=2)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRB_size[cornerCells[i]], pRC_size[cornerCells[i]],candItems);
                        }
                        twoLevel::window::Range_Corners(pR[cornerCells[0]],s, pRC_size[cornerCells[0]], pRD_size[cornerCells[0]],candItems);
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                }
                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);

                                }
                                
                                break;
                        }
                        
                        //////////////////////////////////////////////////////////////////////

                        result = candItems.size();
                        candItems.clear();
                        candItems.reserve(result);

                        tim.start();
                        lsize = insideCells.size();
                        for(int pid = 0; pid < lsize; pid++)
                        {
                            auto cid = insideCells[pid];
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                candItems.push_back(pR[cid][it].id);
                            }
                        }
                        timeInside = tim.stop();

                        tim.start();
                        for(int i = 0; i < 4; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],candItems);
                        }
                        for(int i = 0; i < 2; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],candItems);
                        }
                        for(int i = 0; i < 3; i+=2)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],candItems);
                        }
                        auto cid = cornerCells[0];
                        twoLevel::window::Range_Corners(pR[cid],s, pRC_size[cid], pRD_size[cid],candItems);
                        timeCorner = tim.stop();
                        
                        tim.start();
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                }

                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);

                                }
                                break;
                        
                        }
                        timeBorders = tim.stop();

                        vector<RecordId> resultAfterRefinement;
                        resultAfterRefinement.reserve(result);

                        tim.start();
                        for(auto item : candItems)
                        {
                            if(refinement::refine(s,geo.at(item)))
                                resultAfterRefinement.push_back(item);
                        }
                        double timeRefinement = tim.stop();

                        cout << "twoLevel-Window-Query-Simple-refinement\t" << j << "\t" << candItems.size() << "\t" << resultAfterRefinement.size() << "\t" << candItems.size() << "\t" << timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting << "\t" << 0 << "\t" << timeRefinement << "\t"<<timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting+timeRefinement<< endl;
                    }
                }
                break;
            case REFAVOID:
                tim.start();
                partition::twoLevel::single::PartitionTwoDimensional(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t j = 0; j < size; j++){
                    for(int i = 0; i < NUM_ITERATIONS; i++)
                    {
                        vector<RecordId> candItems;
                        vector<int> insideCells;
                        vector<int> cornerCells;

                        result = 0;
                        auto &s = S[j];
                        
                        tim.start();
                        queryCase = partition::FindRelevantTiles(s,cornerCells,insideCells,runNumPartitionsPerRelation);
                        timeFindingTiles = tim.stop();
                        
                        auto insideCellsSize = insideCells.size();
                        for(int pid = 0; pid < insideCellsSize; pid++)
                        {
                            result += pRA_size[insideCells[pid]];
                        }

                        for(int i = 0; i < 4; i++)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, 0, pRA_size[cornerCells[i]],candItems);
                        }
                        for(int i = 0; i < 2; i++)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRA_size[cornerCells[i]], pRB_size[cornerCells[i]],candItems);
                        }
                        for(int i = 0; i < 3; i+=2)
                        {
                            if(cornerCells[i] != -1)
                                twoLevel::window::Range_Corners(pR[cornerCells[i]],s, pRB_size[cornerCells[i]], pRC_size[cornerCells[i]],candItems);
                        }
                        twoLevel::window::Range_Corners(pR[cornerCells[0]],s, pRC_size[cornerCells[0]], pRD_size[cornerCells[0]],candItems);
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i], candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i], candItems);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i], candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i], candItems);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i], candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i], candItems);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i], candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i], candItems);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i], candItems);
                                }
                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i], candItems);

                                }
                                
                                break;
                        }
                        
                        //////////////////////////////////////////////////////////////////////
                        result = candItems.size();
                        candItems.clear();
                        candItems.reserve(result);

                        tim.start();
                        for(int pid = 0; pid < insideCellsSize; pid++)
                        {
                            auto cid = insideCells[pid];
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                candItems.push_back(pR[cid][it].id);
                            }
                        }
                        timeInside = tim.stop();

                        tim.start();
                        for(int i = 0; i < 4; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],candItems);
                        }
                        for(int i = 0; i < 2; i++)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],candItems);
                        }
                        for(int i = 0; i < 3; i+=2)
                        {
                            auto cid = cornerCells[i];
                            if(cid != -1)
                                twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],candItems);
                        }
                        auto cid = cornerCells[0];
                        twoLevel::window::Range_Corners(pR[cid],s, pRC_size[cid], pRD_size[cid],candItems);
                        timeCorner = tim.stop();
                        
                        tim.start();
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],candItems);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],candItems);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],candItems);
                                }

                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],candItems);

                                }
                                break;
                        
                        }
                        timeBorders = tim.stop();

                        vector<RecordId> refinementObjects;
                        vector<RecordId> resultAfterRefinement;
                        resultAfterRefinement.reserve(result);

                        tim.start();
                        for(auto item : candItems)
                        {
                            if(R[item].xStart >= s.xStart && R[item].xEnd <= s.xEnd)
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            if(R[item].yStart >= s.yStart && R[item].yEnd <= s.yEnd)
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }
                        double timeRefinementTest = tim.stop();

                        tim.start();
                        for(auto item : refinementObjects)
                        {
                            if(refinement::refine(s,geo.at(item)))
                                resultAfterRefinement.push_back(item);
                        }
                        double timeRefinement = tim.stop();
                        cout << "twoLevel-Window-Query-RefAvoid\t" << j << "\t" << candItems.size() << "\t" << resultAfterRefinement.size() << "\t" << refinementObjects.size() << "\t" << timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting << "\t" << timeRefinementTest << "\t" << timeRefinement <<"\t"<<timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting+timeRefinementTest+timeRefinement << endl;
                    }
                }
                break;

            case REFAVOIDPLUS:
                tim.start();
                partition::twoLevel::single::PartitionTwoDimensional(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t j = 0; j < size; j++){
                    for(int i = 0; i < NUM_ITERATIONS; i++)
                    {
                        vector<RecordId> candItems;
                        vector<RecordId> firstCornerCands;
                        vector<RecordId> ABCands;
                        vector<RecordId> ACCands;
                        vector<RecordId> restCands;

                        vector<int> insideCells;
                        vector<int> cornerCells;

                        int restResults = 0;
                        int firstCornerResults = 0;
                        int ABresults = 0;
                        int ACresults=  0;

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
                                restCands.push_back(pR[cid][it].id);
                            }
                        }

                        auto cid = cornerCells[0];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRC_size[cid], pRD_size[cid],firstCornerCands);
                        }

                        cid = cornerCells[1];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],ABCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],ABCands);

                        }

                        cid = cornerCells[2];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],ACCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],ACCands);

                        }
                        
                        cid = cornerCells[3];
                        if(cid != -1)
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],restCands);
                        
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],ABCands);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],ABCands);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],ACCands);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],ACCands);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],ABCands);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],ABCands);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],ACCands);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],ACCands);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],restCands);
                                }

                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],restCands);

                                }
                                break;
                        }
                        
                        result = firstCornerCands.size()+ABCands.size()+ACCands.size()+restCands.size();

                        //////////////////////////////////////////////////////////////////////
                        candItems.reserve(result);

                        firstCornerResults = firstCornerCands.size();
                        firstCornerCands.clear();
                        firstCornerCands.reserve(firstCornerResults);


                        ABresults = ABCands.size();
                        ABCands.clear();
                        ABCands.reserve(ABresults);

                        ACresults = ACCands.size();
                        ACCands.clear();
                        ACCands.reserve(ACresults);

                        restResults = restCands.size();
                        restCands.clear();
                        restCands.reserve(restResults);

                        tim.start();
                        lsize = insideCells.size();
                        for(int pid = 0; pid < lsize; pid++)
                        {
                            auto cid = insideCells[pid];
                            for (auto it = 0; it < pRA_size[cid]; it++) {
                                restCands.push_back(pR[cid][it].id);
                            }
                        }
                        timeInside = tim.stop();

                        tim.start();
                        cid = cornerCells[0];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],firstCornerCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRC_size[cid], pRD_size[cid],firstCornerCands);
                        }

                        cid = cornerCells[1];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],ABCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRA_size[cid], pRB_size[cid],ABCands);

                        }

                        cid = cornerCells[2];
                        if(cid != -1)
                        {
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],ACCands);
                            twoLevel::window::Range_Corners(pR[cid],s, pRB_size[cid], pRC_size[cid],ACCands);

                        }
                        
                        cid = cornerCells[3];
                        if(cid != -1)
                            twoLevel::window::Range_Corners(pR[cid],s, 0, pRA_size[cid],restCands);

                        timeCorner = tim.stop();
                        
                        tim.start();
                        switch (queryCase){
                            case 0: 
                                break;
                            case 1:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],ABCands);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],ABCands);
                                }
                                
                                break;
                            case 2:
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],ACCands);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],ACCands);
                                }
                                
                                break;
                            case 3:
                                
                                for(int i = cornerCells[0]+1; i < cornerCells[1]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],ABCands);
                                    twoLevel::window::Range_B_Class(pR[i],s, pRA_size[i], pRB_size[i],ABCands);
                                }
                                
                                for(int i = cornerCells[0]+runNumPartitionsPerRelation; i < cornerCells[2]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],ACCands);
                                    twoLevel::window::Range_C_Class(pR[i],s, pRB_size[i], pRC_size[i],ACCands);
                                }
                                
                                for( int i = cornerCells[1]+runNumPartitionsPerRelation; i < cornerCells[3]; i+=runNumPartitionsPerRelation)
                                {
                                    twoLevel::window::Range_C_Class(pR[i],s, 0, pRA_size[i],restCands);
                                }

                                for ( int i = cornerCells[2] + 1; i< cornerCells[3]; i++)
                                {
                                    twoLevel::window::Range_B_Class(pR[i],s, 0, pRA_size[i],restCands);

                                }
                                break;
                        
                        }
                        
                        timeBorders = tim.stop();

                        vector<RecordId> resultAfterRefinement;
                        resultAfterRefinement.reserve(result);
                        
                        vector<RecordId> refinementObjects;

                        tim.start();
                        for(auto item : firstCornerCands)
                        {
                            if(R[item].xStart >= s.xStart && R[item].xEnd <= s.xEnd)
                            {   
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            if(R[item].yStart >= s.yStart && R[item].yEnd <= s.yEnd)
                            {   
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }

                        for(auto item : ABCands)
                        {
                            if(R[item].xEnd <= s.xEnd || (R[item].yStart >= s.yStart && R[item].yEnd <= s.yEnd))
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }

                        for(auto item : ACCands)
                        {
                            if(R[item].yEnd <= s.yEnd || (R[item].xStart >= s.xStart && R[item].xEnd <= s.xEnd))
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }

                        for(auto item : restCands)
                        {
                            if((R[item].xEnd <= s.xEnd) || (R[item].yEnd <= s.yEnd))
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }
                        double timeRefinementTest = tim.stop();

                        tim.start();
                        for(auto item : refinementObjects)
                        {
                            if(refinement::refine(s,geo.at(item)))
                                resultAfterRefinement.push_back(item);

                        }
                        double timeRefinement = tim.stop();
                        unsigned long long newRes = firstCornerCands.size() + ABCands.size() + ACCands.size() + restCands.size();
                        cout << "twoLevel-Window-Query-RefAvoid+\t" << j << "\t" << newRes << "\t" << resultAfterRefinement.size() << "\t" << refinementObjects.size() << "\t" << timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting << "\t" << timeRefinementTest << "\t" << timeRefinement<<"\t"<<timeFindingTiles+timeInside+timeCorner+timeBorders+timeCounting+timeRefinementTest+timeRefinement << endl;
                    }
                }
                break;
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
        refinement::normalizeGeometries(minX,minY,maxExtend,geo);


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


        switch(refinementCase)
        {  
            case SIMPLE_REFINEMENT:
                tim.start();
                partition::twoLevel::single::PartitionTwoDimensional(R, pR,pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                timeIndexingOrPartitioning = tim.stop();

                for (size_t j = 0; j < S.numRecords; j++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) 
                    {
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

                        vector<RecordId> resultAfterRefinement;
                        resultAfterRefinement.reserve(resultItemsIds.size());

                        Coord centerX = (s.xStart + s.xEnd)/2;
                        Coord centerY = (s.yStart + s.yEnd)/2;
                        Coordinates centerPoint(centerX,centerY);

                        tim.start();
                        for(auto item : resultItemsIds)
                        {
                            if(refinement::refineDisc(centerPoint,geo.at(item),epsilon))
                                resultAfterRefinement.push_back(item);
                        }
                        double timeRefinement = tim.stop();

                        cout << "twoLevel-Disk-Query-Simple-refinement\t" << j << "\t" << newRes << "\t" << resultAfterRefinement.size() << "\t" << newRes << "\t" << timeInside+timeCorner+timeBorders+timeCounting+timeFindingTiles << "\t" << 0 << "\t" << timeRefinement<<"\t"<<timeInside+timeCorner+timeBorders+timeCounting+timeFindingTiles+timeRefinement << endl;
                    
                    }
                }
                break;

            case REFAVOID:
                tim.start();
                partition::twoLevel::single::PartitionTwoDimensional(R, pR,pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);
                timeIndexingOrPartitioning = tim.stop();

                
                for (size_t j = 0; j < S.numRecords; j++){
                    for(int o = 0; o < NUM_ITERATIONS; o++) 
                    {
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
                        
                        vector<RecordId> refinementObjects;

                        Coord centerX = (s.xStart + s.xEnd)/2;
                        Coord centerY = (s.yStart + s.yEnd)/2;
                        Coordinates centerPoint(centerX,centerY);

                        tim.start();

                        vector<RecordId> resultAfterRefinement;
                        resultAfterRefinement.reserve(resultItemsIds.size());

                        for(auto item : resultItemsIds)
                        {
                            int counter=0;
                            if(twoLevel::disk::distPointToMbr(R[item].xStart,R[item].yStart,centerPoint) <= epsilon)
                                counter++;
                            if(twoLevel::disk::distPointToMbr(R[item].xEnd,R[item].yStart,centerPoint) <= epsilon)
                                counter++;
                            if(twoLevel::disk::distPointToMbr(R[item].xEnd,R[item].yEnd,centerPoint) <= epsilon)
                                counter++;
                            if(twoLevel::disk::distPointToMbr(R[item].xStart,R[item].yEnd,centerPoint) <= epsilon)
                                counter++;

                            if(counter > 1)
                            {
                                resultAfterRefinement.push_back(item);
                                continue;
                            }
                            refinementObjects.push_back(item);
                        }
                        double timeRefinementTest = tim.stop();

                    

                        tim.start();
                        for(auto item : refinementObjects)
                        {
                            if(refinement::refineDisc(centerPoint,geo.at(item),epsilon))
                                resultAfterRefinement.push_back(item);
                        }
                        double timeRefinement = tim.stop();

                        cout << "twoLevel-Disk-Query-RefAvoid\t" << j << "\t" << newRes << "\t" << resultAfterRefinement.size() << "\t" << (double)refinementObjects.size() << "\t" << timeInside+timeCorner+timeBorders+timeCounting+timeFindingTiles << "\t" << timeRefinementTest << "\t" << timeRefinement<<"\t"<<timeInside+timeCorner+timeBorders+timeCounting+timeFindingTiles+timeRefinementTest+timeRefinement<< endl;
                    }
                }
                break;
        }
        break;
    }
    
}
