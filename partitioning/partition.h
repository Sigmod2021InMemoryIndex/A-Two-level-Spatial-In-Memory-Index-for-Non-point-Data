
namespace partition
{
    int myQuotient(double numer, double denom) {
        return int(numer/denom + EPS);
    };


    double myRemainder(double numer, double denom, int q) {
        double rem = double(numer - q*denom);

        return ((abs(rem) < EPS) ? 0: rem);
    };


    int getCellId(int x, int y, int numCellsPerDimension) {
        return (y * numCellsPerDimension + x);
    };


    int findReferenceCell(double x, double y, double cellExtent, int numCellsPerDimension) {
        int xInt,yInt;

        xInt = (x + EPS)/cellExtent;
        yInt = (y + EPS)/cellExtent;

        return (yInt * numCellsPerDimension + xInt);
    };

    int FindRelevantTiles(Record &s,vector<int> &cornerCells,vector<int> &insideCells,int runNumPartitionsPerRelation)
    {
        double xStartCell, yStartCell, xEndCell, yEndCell;

        int queryCase;
        int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
        Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
        xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
        yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
        auto xEnd = myRemainder(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
        auto yEnd = myRemainder(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));
        

        if (s.xEnd + EPS >= 1) {
            xEndCell = runNumPartitionsPerRelation - 1;
        }
        else {
            xEndCell = myQuotient(s.xEnd + EPS, partitionExtent);
        }

        if (s.yEnd + EPS >= 1) {
            yEndCell = runNumPartitionsPerRelation - 1;
        }
        else {
            yEndCell = myQuotient(s.yEnd + EPS, partitionExtent);
        }

        int bottomLeft = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
        int bottomRight = getCellId(xEndCell, yStartCell, runNumPartitionsPerRelation);
        int topLeft = getCellId(xStartCell, yEndCell, runNumPartitionsPerRelation);
        int topRight = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);
        
        
        if (bottomLeft == topLeft && bottomLeft == bottomRight ){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            queryCase = 0 ;
        }
        else if (bottomLeft == topLeft){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(bottomRight);
            cornerCells.push_back(-1);
            cornerCells.push_back(-1);
            queryCase = 1;
        }
        else if ( bottomLeft == bottomRight){
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(-1);
            cornerCells.push_back(topLeft);
            cornerCells.push_back(-1);
            queryCase = 2;
        }
        else{
        
            cornerCells.push_back(bottomLeft);
            cornerCells.push_back(bottomRight);
            cornerCells.push_back(topLeft);
            cornerCells.push_back(topRight);
            queryCase = 3;
        }

        int rightLimit = bottomRight;
        int start = bottomLeft;

        int ii = start; 
        while( ii<=topRight)
        {
            if(ii <= rightLimit)
            {
                if(ii != start && ii != rightLimit && start != bottomLeft && rightLimit != topRight)  
                {
                    insideCells.push_back(ii);
                }
                ii++;
            }
            else{
                rightLimit +=runNumPartitionsPerRelation;
                start += runNumPartitionsPerRelation;
                ii=start;
            }
        }
        
        return queryCase;
    }


    namespace oneLevel
    {
        namespace single
        {
            void PartitionUniform(const Relation& R, Relation *pR, int runNumPartitionsPerRelation)
            {
                int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
                Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
                size_t *pR_size = new size_t[runNumPartitions];

                memset(pR_size, 0, runNumPartitions*sizeof(size_t));
   
                double xStartCell, yStartCell, xEndCell, yEndCell;
                int firstCell, lastCell;
   
                for (size_t i = 0; i < R.numRecords; i++){
                    auto &r = R[i];
                    // Determine cell for (rec.xStart, rec.yStart)
                    xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
                    yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }


                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation-1;
                    }
                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }

                    lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

                    // Put record in cells.
                    if (firstCell == lastCell) {
                        pR_size[firstCell]++;

                    }
                    else {
                        pR_size[firstCell]++;

                        int cellNumber;
                        for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                            if ( i != xStartCell){
                                cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                                pR_size[cellNumber]++;
                            }
                            for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                cellNumber = getCellId(i, j, runNumPartitionsPerRelation);

                                pR_size[cellNumber]++;
                            }
                        }
                    }
                }

                for (int i = 0; i < runNumPartitions; i++){
                    pR[i].reserve(pR_size[i]);
                    pR[i].numRecords = pR_size[i];
                }

                for (size_t i = 0; i < R.numRecords; i++){
                    
                    auto &r = R[i];

                    xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
                    yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }

                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }

                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation-1;
                    }

                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }

                    lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

                    // Put record in cells.
                    if (firstCell == lastCell) {
                        pR[firstCell].push_back(r);

                    }
                    else {

                        pR[firstCell].push_back(r);

                        int cellNumber;
                        for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                            if ( i != xStartCell){
                                cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                                pR[cellNumber].push_back(r);
                            }
                            for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                                pR[cellNumber].push_back(r);
                            }
                        }
                    }        
                }
            }

            void PartitionTwoDimensional(const Relation& R, Relation *pR, int runNumPartitionsPerRelation)
            {
                int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
                PartitionUniform(R, pR, runNumPartitionsPerRelation);
            };
        } 
        namespace parallel
        {

            
        } 
        
    } 

    namespace twoLevel
    {
        namespace single
        {

            void PartitionUniform(const Relation& R, Relation *pR, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation)
            {
                int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
                Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

                double xStartCell, yStartCell, xEndCell, yEndCell;
                int firstCell, lastCell;
                Timer tim;
                double timepR = 0, timeDecomp = 0;

                for (size_t i = 0; i < R.numRecords; i++){
                    auto &r = R[i];

                    // Determine cell for (rec.xStart, rec.yStart)
                    xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                    yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }

                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }

                    lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                    int x = 0, y = 0;

                    // Put record in cells.
                    if (firstCell == lastCell) {
                        pRA_size[firstCell]++;
                    }
                    else {
                        pRA_size[firstCell]++;
                        int cellNumber;
                        for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                            if ( i != xStartCell){
                                cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                                pRC_size[cellNumber]++;
                            }
                            for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                                if( i == xStartCell){
                                    pRB_size[cellNumber]++;
                                }
                                else{
                                    pRD_size[cellNumber]++;
                                }
                            }
                        }
                    }
                }

                int counter = 0;
                for (int i = 0; i < runNumPartitions; i++){
                    counter = pRA_size[i] + pRB_size[i] + pRC_size[i] + pRD_size[i] ;
                    pR[i].resize(counter);
                    pR[i].numRecords = counter;
                }

                for (int i = 0; i < runNumPartitions; i++){
                    pRD_size[i] = pRC_size[i] + pRB_size[i] + pRA_size[i];
                    pRC_size[i] = pRB_size[i] + pRA_size[i];
                    pRB_size[i] = pRA_size[i];
                    pRA_size[i] = 0;           
                    
                }

                for (size_t i = 0; i < R.numRecords; i++){
                    auto &r = R[i];

                    xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                    yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                    firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                    // Determine cell for (rec.xEnd, rec.yEnd)
                    auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                    auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                    if (r.xEnd + EPS >= 1) {
                        xEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                    }

                    if (r.yEnd + EPS >= 1) {
                        yEndCell = runNumPartitionsPerRelation - 1;
                    }
                    else {
                        yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                    }
                    lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                    int x = 0 , y = 0;

                    // Put record in cells.
                    if (firstCell == lastCell) {

                        pR[firstCell][pRA_size[firstCell]] = r;
                        pRA_size[firstCell] = pRA_size[firstCell] + 1;
                        
                    }
                    else {

                        pR[firstCell][pRA_size[firstCell]] = r;
                        pRA_size[firstCell] = pRA_size[firstCell] + 1;

                        int cellNumber;
                        for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                            if ( i != xStartCell){
                                cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                                pR[cellNumber][pRC_size[cellNumber]] = r;
                                pRC_size[cellNumber] = pRC_size[cellNumber] + 1;

                            }
                            for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                                if( i == xStartCell){

                                    pR[cellNumber][pRB_size[cellNumber]] = r;
                                    pRB_size[cellNumber] = pRB_size[cellNumber] + 1 ;

                                }
                                else{

                                    pR[cellNumber][pRD_size[cellNumber]] = r;
                                    pRD_size[cellNumber] = pRD_size[cellNumber] + 1 ;
                                }
                            }
                        }
                    }
                }
            };


            void PartitionTwoDimensional(Relation& R, Relation *pR, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation)
            {
                int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
                PartitionUniform(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);            
            };
            
        } 
        namespace parallel
        {

            void Partition_Range(const Relation& R, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation, int insertPercentage, int runNumThreads)
            {
                int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
                Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
                size_t **tmp_pRA_size, **tmp_pRB_size,  **tmp_pRC_size, **tmp_pRD_size;
                size_t chunkSizeR = R.numRecords/runNumThreads;
                
                int extend_vec = ceil(insertPercentage*R.numRecords);
                #pragma omp parallel // shared(pR, pS)
                {
                    double xStartCell, yStartCell, xEndCell, yEndCell;
                    int firstCell, lastCell;
                    
                    #pragma omp single
                    {
                        runNumThreads = omp_get_num_threads();
                        tmp_pRA_size = new size_t*[runNumThreads];
                        tmp_pRB_size = new size_t*[runNumThreads];
                        tmp_pRC_size = new size_t*[runNumThreads];
                        tmp_pRD_size = new size_t*[runNumThreads];
                    }

                    #pragma omp for
                    for (int th = 0; th < runNumThreads; th++)
                    {
                        tmp_pRA_size[th] = new size_t[runNumPartitions];
                        tmp_pRB_size[th] = new size_t[runNumPartitions];
                        tmp_pRC_size[th] = new size_t[runNumPartitions];
                        tmp_pRD_size[th] = new size_t[runNumPartitions];
                        
                        memset(tmp_pRA_size[th], 0, runNumPartitions*sizeof(size_t));
                        memset(tmp_pRB_size[th], 0, runNumPartitions*sizeof(size_t));
                        memset(tmp_pRC_size[th], 0, runNumPartitions*sizeof(size_t));
                        memset(tmp_pRD_size[th], 0, runNumPartitions*sizeof(size_t));
                    }
                    
                    #pragma omp for reduction(+ : pRA_size[:runNumPartitions], pRB_size[:runNumPartitions], pRC_size[:runNumPartitions], pRD_size[:runNumPartitions]) schedule(static, chunkSizeR)
                    for (size_t i = 0; i < R.numRecords; i++){
                        auto &r = R[i];
                        auto th = omp_get_thread_num();

                        // Determine cell for (rec.xStart, rec.yStart)
                        xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                        yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                        firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                        // Determine cell for (rec.xEnd, rec.yEnd)
                        auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                        auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));


                        if (r.xEnd + EPS >= 1) {
                            xEndCell = runNumPartitionsPerRelation - 1;
                        }
                        else {
                            xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                        }

                        if (r.yEnd + EPS >= 1) {
                            yEndCell = runNumPartitionsPerRelation - 1;
                        }
                        else {
                            yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                        }

                        lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                        int x = 0, y = 0;

                    
                        // Put record in cells.
                        if (firstCell == lastCell) {
                            pRA_size[firstCell]++;
                            tmp_pRA_size[th][firstCell]++;
                        }
                        else {
                            pRA_size[firstCell]++;
                            tmp_pRA_size[th][firstCell]++;
                            
                            int cellNumber;
                            for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                                if ( i != xStartCell){
                                    cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                                    
                                    pRC_size[cellNumber]++;
                                    tmp_pRC_size[th][cellNumber]++;
                                }
                                for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                    cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                                    if( i == xStartCell){
                                        
                                        pRB_size[cellNumber]++;
                                        tmp_pRB_size[th][cellNumber]++;
                                    }
                                    else{
                                        pRD_size[cellNumber]++;
                                        tmp_pRD_size[th][cellNumber]++;
                                    }
                                }
                            }
                        }
                    }

                    #pragma omp for 
                    for (int pid = 0; pid < runNumPartitions; pid++)
                    {
                        pRA[pid].resize(pRA_size[pid] + extend_vec);
                        pRB[pid].resize(pRB_size[pid] + extend_vec);
                        pRC[pid].resize(pRC_size[pid] + extend_vec);
                        pRD[pid].resize(pRD_size[pid] + extend_vec);
                    }

                    #pragma omp single
                    for (int pid = 0; pid < runNumPartitions; pid++)
                    {
                        auto prevr = 0;
                        for (int th = 0; th < runNumThreads; th++)
                        {
                            auto tmpr = tmp_pRA_size[th][pid];
                            tmp_pRA_size[th][pid] = prevr;
                            prevr += tmpr;
                        }
                        
                        auto prevr2 = 0;
                        for ( int th = 0 ; th < runNumThreads; th++){
                            
                            auto tmpr2 = tmp_pRB_size[th][pid];
                            tmp_pRB_size[th][pid] = prevr2;
                            prevr2 += tmpr2;
                        }
                        
                        auto prevr3 = 0;
                        for ( int th = 0 ; th < runNumThreads; th++){
                            
                            auto tmpr3 = tmp_pRC_size[th][pid];
                            tmp_pRC_size[th][pid] = prevr3;
                            prevr3 += tmpr3;
                        }
                    
                        auto prevr4 = 0;
                        for ( int th = 0 ; th < runNumThreads; th++){
                            
                            auto tmpr4 = tmp_pRD_size[th][pid];
                            tmp_pRD_size[th][pid] = prevr4;
                            prevr4 += tmpr4;
                            
                        }
                    }

                    #pragma omp for schedule(static, chunkSizeR)                      
                    for (size_t i = 0; i < R.numRecords; i++)
                    {
                        auto &r = R[i];
                        auto th = omp_get_thread_num();
                        
                        xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                        yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                        firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                        // Determine cell for (rec.xEnd, rec.yEnd)
                        auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                        auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                        if (r.xEnd + EPS >= 1) {
                            xEndCell = runNumPartitionsPerRelation - 1;
                        }
                        else {
                            xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                        }

                        if (r.yEnd + EPS >= 1) {
                            yEndCell = runNumPartitionsPerRelation - 1;
                        }
                        else {
                            yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                        }
                        lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                        int x = 0 , y = 0;

                        // Put record in cells.
                        if (firstCell == lastCell) {
                            pRA[firstCell][tmp_pRA_size[th][firstCell]] = r;
                            tmp_pRA_size[th][firstCell]++;
                        }
                        else {
                            
                            pRA[firstCell][tmp_pRA_size[th][firstCell]] = r;
                            
                            
                            tmp_pRA_size[th][firstCell]++;

                            int cellNumber;
                            for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                                if ( i != xStartCell){
                                    cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                                    pRC[cellNumber][tmp_pRC_size[th][cellNumber]] = r; 
                                    tmp_pRC_size[th][cellNumber]++;                                       
                                }
                                for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                                    cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                                    if( i == xStartCell){
                                        pRB[cellNumber][tmp_pRB_size[th][cellNumber]] = r;
                                        tmp_pRB_size[th][cellNumber]++;
                                    }
                                    else{
                                        pRD[cellNumber][tmp_pRD_size[th][cellNumber]] = r;
                                        tmp_pRD_size[th][cellNumber]++;
                                    }
                                }
                            }
                        }
                    }
                    
                    #pragma omp for
                    for (int th = 0; th < runNumThreads; th++)
                    {
                        delete[] tmp_pRA_size[th];
                        delete[] tmp_pRB_size[th];
                        delete[] tmp_pRC_size[th];
                        delete[] tmp_pRD_size[th];
                    }
                }

                delete[] tmp_pRA_size;
                delete[] tmp_pRB_size;
                delete[] tmp_pRC_size;
                delete[] tmp_pRD_size;
                    
            };
                
            void PartitionTwoDimensionalRange(Relation& R, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation, int insertPercentage, int runNumThreads)
            {
                int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
                Partition_Range(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation, insertPercentage, runNumThreads);                        
            };
        } 
    } 

    namespace twoLevelPlus
    {

        void PartitionUniform(const Relation& R, Relation *pR, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, vector<Decompose> *pRXStart, vector<Decompose> *pRXEnd, vector<Decompose> *pRYStart, vector<Decompose> *pRYEnd, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;
            Timer tim;
            double timepR = 0, timeDecomp = 0;

            for (size_t i = 0; i < R.numRecords; i++){
                auto &r = R[i];

                // Determine cell for (rec.xStart, rec.yStart)
                xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }

                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }

                lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                int x = 0, y = 0;

                // Put record in cells.
                if (firstCell == lastCell) {
                    pRA_size[firstCell]++;
                }
                else {
                    pRA_size[firstCell]++;
                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                            pRC_size[cellNumber]++;
                        }
                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                            if( i == xStartCell){
                                pRB_size[cellNumber]++;
                            }
                            else{
                                pRD_size[cellNumber]++;
                            }
                        }
                    }
                }
            }

            int counter = 0;
            for (int i = 0; i < runNumPartitions; i++){
                counter = pRA_size[i] + pRB_size[i] + pRC_size[i] + pRD_size[i] ;
                pR[i].resize(counter);
                pR[i].numRecords = counter;
                
                pRXStart[i].resize(counter);
                pRXEnd[i].resize(counter);
                pRYStart[i].resize(counter);
                pRYEnd[i].resize(counter);
                
            }

            for (int i = 0; i < runNumPartitions; i++){
                pRD_size[i] = pRC_size[i] + pRB_size[i] + pRA_size[i];
                pRC_size[i] = pRB_size[i] + pRA_size[i];
                pRB_size[i] = pRA_size[i];
                pRA_size[i] = 0;           
                
            }

            for (size_t i = 0; i < R.numRecords; i++){
                auto &r = R[i];

                xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }

                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }
                lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                int x = 0 , y = 0;

                // Put record in cells.
                if (firstCell == lastCell) {

                    pR[firstCell][pRA_size[firstCell]] = r;
                    
                    pRXStart[firstCell][pRA_size[firstCell]].id = r.id;
                    pRXStart[firstCell][pRA_size[firstCell]].value = r.xStart;
                    
                    pRXEnd[firstCell][pRA_size[firstCell]].id = r.id;
                    pRXEnd[firstCell][pRA_size[firstCell]].value = r.xEnd;                            
                    
                    pRYStart[firstCell][pRA_size[firstCell]].id = r.id;
                    pRYStart[firstCell][pRA_size[firstCell]].value = r.yStart;                            
                    
                    pRYEnd[firstCell][pRA_size[firstCell]].id = r.id;
                    pRYEnd[firstCell][pRA_size[firstCell]].value = r.yEnd;
                    
                    pRA_size[firstCell] = pRA_size[firstCell] + 1;
                    
                }
                else {

                    pR[firstCell][pRA_size[firstCell]] = r;
                    
                    pRXStart[firstCell][pRA_size[firstCell]].id = r.id;
                    pRXStart[firstCell][pRA_size[firstCell]].value = r.xStart;
                    
                    pRXEnd[firstCell][pRA_size[firstCell]].id = r.id;
                    pRXEnd[firstCell][pRA_size[firstCell]].value = r.xEnd;                            
                    
                    pRYStart[firstCell][pRA_size[firstCell]].id = r.id;
                    pRYStart[firstCell][pRA_size[firstCell]].value = r.yStart;                            
                    
                    pRYEnd[firstCell][pRA_size[firstCell]].id = r.id;
                    pRYEnd[firstCell][pRA_size[firstCell]].value = r.yEnd;
                    
                    pRA_size[firstCell] = pRA_size[firstCell] + 1;

                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                            pR[cellNumber][pRC_size[cellNumber]] = r;
                            
                            pRXStart[cellNumber][pRC_size[cellNumber]].id = r.id;
                            pRXStart[cellNumber][pRC_size[cellNumber]].value = r.xStart;

                            pRXEnd[cellNumber][pRC_size[cellNumber]].id = r.id;
                            pRXEnd[cellNumber][pRC_size[cellNumber]].value = r.xEnd;                            

                            pRYStart[cellNumber][pRC_size[cellNumber]].id = r.id;
                            pRYStart[cellNumber][pRC_size[cellNumber]].value = r.yStart;                            

                            pRYEnd[cellNumber][pRC_size[cellNumber]].id = r.id;
                            pRYEnd[cellNumber][pRC_size[cellNumber]].value = r.yEnd;
                            
                            pRC_size[cellNumber] = pRC_size[cellNumber] + 1;

                        }
                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                            if( i == xStartCell){

                                pR[cellNumber][pRB_size[cellNumber]] = r;
                                
                                pRXStart[cellNumber][pRB_size[cellNumber]].id = r.id;
                                pRXStart[cellNumber][pRB_size[cellNumber]].value = r.xStart;

                                pRXEnd[cellNumber][pRB_size[cellNumber]].id = r.id;
                                pRXEnd[cellNumber][pRB_size[cellNumber]].value = r.xEnd;                            

                                pRYStart[cellNumber][pRB_size[cellNumber]].id = r.id;
                                pRYStart[cellNumber][pRB_size[cellNumber]].value = r.yStart;                            

                                pRYEnd[cellNumber][pRB_size[cellNumber]].id = r.id;
                                pRYEnd[cellNumber][pRB_size[cellNumber]].value = r.yEnd;
                                
                                pRB_size[cellNumber] = pRB_size[cellNumber] + 1 ;

                            }
                            else{

                                pR[cellNumber][pRD_size[cellNumber]] = r;
                                
                                pRXStart[cellNumber][pRD_size[cellNumber]].id = r.id;
                                pRXStart[cellNumber][pRD_size[cellNumber]].value = r.xStart;

                                pRXEnd[cellNumber][pRD_size[cellNumber]].id = r.id;
                                pRXEnd[cellNumber][pRD_size[cellNumber]].value = r.xEnd;                            

                                pRYStart[cellNumber][pRD_size[cellNumber]].id = r.id;
                                pRYStart[cellNumber][pRD_size[cellNumber]].value = r.yStart;                            

                                pRYEnd[cellNumber][pRD_size[cellNumber]].id = r.id;
                                pRYEnd[cellNumber][pRD_size[cellNumber]].value = r.yEnd;
                                
                                pRD_size[cellNumber] = pRD_size[cellNumber] + 1 ;
                            }
                        }
                    }
                }
            }

        };

        void PartitionTwoDimensional(Relation& R, Relation *pR, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, vector<Decompose> *pRXStart, vector<Decompose> *pRXEnd, vector<Decompose> *pRYStart, vector<Decompose> *pRYEnd, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
            PartitionUniform(R, pR, pRA_size, pRB_size, pRC_size, pRD_size, pRXStart, pRXEnd, pRYStart, pRYEnd, runNumPartitionsPerRelation);            
        };
    }

    namespace batch
    {

        void FindRelevantTiles(Relation& S, Relation *pSIns, Relation *pSCorDR, Relation *pSCorDL, Relation *pSCorUR, Relation *pSCorUL, Relation *pSBorX, Relation *pSBorY, Relation *pSBorB, Relation *pSBorC, size_t *pSIns_size, size_t *pSCorDR_size, size_t *pSCorDL_size, size_t *pSCorUR_size, size_t *pSCorUL_size, size_t *pSBorX_size, size_t *pSBorY_size, size_t *pSBorB_size, size_t *pSBorC_size, bool * checkTile, int runNumPartitionsPerRelation){
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;
            int queryCase;

            int querySize = S.numRecords;

            for (size_t i = 0; i < querySize; i++){ 
                auto &s = S[i];
                xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
                auto xEnd = myRemainder(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));

                if (s.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(s.xEnd + EPS, partitionExtent);
                }

                if (s.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(s.yEnd + EPS, partitionExtent);
                }

                int downLeft = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
                int downRight = getCellId(xEndCell, yStartCell, runNumPartitionsPerRelation);
                int upLeft = getCellId(xStartCell, yEndCell, runNumPartitionsPerRelation);
                int upRight = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                if (downLeft == upLeft && downLeft == downRight ){    
                    pSCorDL_size[downLeft]++;

                    pSCorDL_size[downLeft]++;
                    checkTile[downLeft] = true;
                }
                else if (downLeft == upLeft){
                    pSCorDL_size[downLeft]++;
                    checkTile[downLeft] = true;


                    pSCorDR_size[downRight]++;
                    checkTile[downRight] = true;

                    for(int j = downLeft + 1; j < downRight; j++)
                    {
                        pSBorX_size[j]++;
                        checkTile[j] = true;
                    } 

                }
                else if ( downLeft == downRight){

                    pSCorDL_size[downLeft]++;
                    checkTile[downLeft] = true;                                

                    pSCorUL_size[upLeft]++;
                    checkTile[upLeft] = true;

                    for(int j = downLeft+runNumPartitionsPerRelation; j < upLeft; j += runNumPartitionsPerRelation)
                    {
                        pSBorY_size[j]++;
                        checkTile[j] = true;
                    }

                }
                else{
                    pSCorDL_size[downLeft]++;
                    checkTile[downLeft] = true;

                    pSCorDR_size[downRight]++;
                    checkTile[downRight] = true;

                    pSCorUL_size[upLeft]++;
                    checkTile[upLeft] = true;


                    pSCorUR_size[upRight]++;
                    checkTile[upRight] = true;;

                    for(int j = downLeft+1; j < downRight; j++)
                    {
                        pSBorX_size[j]++;
                        checkTile[j] = true;
                    }

                    for(int j = downLeft+runNumPartitionsPerRelation; j < upLeft; j+=runNumPartitionsPerRelation)
                    {
                        pSBorY_size[j]++;
                        checkTile[j] = true;
                    }

                    for( int j = downRight+runNumPartitionsPerRelation; j < upRight; j+=runNumPartitionsPerRelation)
                    {
                        pSBorC_size[j]++;
                        checkTile[j] = true;
                    }
                    for ( int j = upLeft + 1; j< upRight; j++)
                    {
                        pSBorB_size[j]++;
                        checkTile[j] = true;
                    }
                }

                int rightLimit = downRight;
                int start = downLeft;

                int ii = start; 
                while( ii<=upRight)
                {
                    if(ii <= rightLimit)
                    {
                        if(ii != start && ii != rightLimit && start != downLeft && rightLimit != upRight)  
                        {
                            pSIns_size[ii]++;
                            checkTile[ii] = true;
                        }
                        ii++;
                    }
                    else{
                        rightLimit +=runNumPartitionsPerRelation;
                        start += runNumPartitionsPerRelation;
                        ii=start;
                    }
                }
            }

            for (int pid = 0; pid < runNumPartitions; pid++)
            {

                pSIns[pid].reserve(pSIns_size[pid]);
                pSCorDR[pid].reserve(pSCorDR_size[pid]);
                pSCorDL[pid].reserve(pSCorDL_size[pid]);
                pSCorUR[pid].reserve(pSCorUR_size[pid]);
                pSCorUL[pid].reserve(pSCorUL_size[pid]);
                pSBorX[pid].reserve(pSBorX_size[pid]);
                pSBorY[pid].reserve(pSBorY_size[pid]);
                pSBorB[pid].reserve(pSBorB_size[pid]);
                pSBorC[pid].reserve(pSBorC_size[pid]);
            }

            for (size_t i = 0; i < querySize; i++){ 
                auto &s = S[i];
                xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
                auto xEnd = myRemainder(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));


                if (s.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(s.xEnd + EPS, partitionExtent);
                }

                if (s.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(s.yEnd + EPS, partitionExtent);
                }

                int downLeft = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
                int downRight = getCellId(xEndCell, yStartCell, runNumPartitionsPerRelation);
                int upLeft = getCellId(xStartCell, yEndCell, runNumPartitionsPerRelation);
                int upRight = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                if (downLeft == upLeft && downLeft == downRight ){    
                    pSCorDL[downLeft].push_back(s);
                }
                else if (downLeft == upLeft){
                    pSCorDL[downLeft].push_back(s);                                

                    pSCorDR[downRight].push_back(s);                                

                    for(int j = downLeft + 1; j < downRight; j++)
                    {
                        pSBorX[j].push_back(s);
                    } 

                }
                else if ( downLeft == downRight){

                    pSCorDL[downLeft].push_back(s);                                

                    pSCorUL[upLeft].push_back(s);                                

                    for(int j = downLeft+runNumPartitionsPerRelation; j < upLeft; j += runNumPartitionsPerRelation)
                    {
                        pSBorY[j].push_back(s);                                      
                    }

                }
                else{
                    pSCorDL[downLeft].push_back(s);

                    pSCorDR[downRight].push_back(s);

                    pSCorUL[upLeft].push_back(s);                                

                    pSCorUR[upRight].push_back(s);                                

                    for(int j = downLeft+1; j < downRight; j++)
                    {
                        pSBorX[j].push_back(s);
                    }

                    for(int j = downLeft+runNumPartitionsPerRelation; j < upLeft; j+=runNumPartitionsPerRelation)
                    {
                        pSBorY[j].push_back(s);
                    }

                    for( int j = downRight+runNumPartitionsPerRelation; j < upRight; j+=runNumPartitionsPerRelation)
                    {
                        pSBorC[j].push_back(s);
                    }
                    for ( int j = upLeft + 1; j< upRight; j++)
                    {
                        pSBorB[j].push_back(s);
                    }
                }

                int rightLimit = downRight;
                int start = downLeft;

                int ii = start; 
                while( ii<=upRight)
                {
                    if(ii <= rightLimit)
                    {
                        if(ii != start && ii != rightLimit && start != downLeft && rightLimit != upRight)  
                        {
                            pSIns[ii].push_back(s);
                        }
                        ii++;
                    }
                    else{
                        rightLimit +=runNumPartitionsPerRelation;
                        start += runNumPartitionsPerRelation;
                        ii=start;
                    }
                }
            }                  
        }
    
    }

    namespace update
    {
        void UpdateCounters(Record &r, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation){
            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            
            xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
            yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
            firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

            // Determine cell for (rec.xEnd, rec.yEnd)
            auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

            if (r.xEnd + EPS >= 1) {
                xEndCell = runNumPartitionsPerRelation - 1;
            }
            else {
                xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
            }

            if (r.yEnd + EPS >= 1) {
                yEndCell = runNumPartitionsPerRelation - 1;
            }
            else {
                yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
            }
            lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

            int x = 0 , y = 0;

            // Put record in cells.
            if (firstCell == lastCell) {
                pRA_size[firstCell] = pRA_size[firstCell] + 1;
            }
            else {
                pRA_size[firstCell] = pRA_size[firstCell] + 1;

                int cellNumber;
                for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                    if ( i != xStartCell){
                        cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                        pRC_size[cellNumber] = pRC_size[cellNumber] + 1;
                    }
                    for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                        cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                        if( i == xStartCell){
                            pRB_size[cellNumber] = pRB_size[cellNumber] + 1 ;
                        }
                        else{
                            pRD_size[cellNumber] = pRD_size[cellNumber] + 1 ;
                        }
                    }
                }
            }
        }

        void UpdateCounters( Record & r, Relation *pR, size_t *pR_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;

            xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
            yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
            firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

            // Determine cell for (rec.xEnd, rec.yEnd)
            auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

            if (r.xEnd + EPS >= 1) {
                xEndCell = runNumPartitionsPerRelation - 1;
            }

            else {
                xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
            }

            if (r.yEnd + EPS >= 1) {
                yEndCell = runNumPartitionsPerRelation-1;
            }

            else {
                yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
            }

            lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

            // Put record in cells.
            if (firstCell == lastCell) {
                pR_size[firstCell] = pR_size[firstCell] + 1;              
            }
            else {
                pR_size[firstCell] = pR_size[firstCell] + 1;

                int cellNumber;
                for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                    if ( i != xStartCell){
                        cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                        pR_size[cellNumber] = pR_size[cellNumber] + 1;
                    }
                    for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                        cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                        pR_size[cellNumber] = pR_size[cellNumber] + 1;
                    }
                }
            }        
            
        }
    
        void insert( Record & r, Relation *pR, size_t *pR_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;

            xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
            yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
            firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

            // Determine cell for (rec.xEnd, rec.yEnd)
            auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

            if (r.xEnd + EPS >= 1) {
                xEndCell = runNumPartitionsPerRelation - 1;
            }

            else {
                xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
            }

            if (r.yEnd + EPS >= 1) {
                yEndCell = runNumPartitionsPerRelation-1;
            }

            else {
                yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
            }

            lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

            // Put record in cells.
            if (firstCell == lastCell) {
                pR[firstCell][pR_size[firstCell]] = r;
                pR_size[firstCell] = pR_size[firstCell] + 1;         
            }
            else {
                pR[firstCell][pR_size[firstCell]] = r;
                pR_size[firstCell] = pR_size[firstCell] + 1;

                int cellNumber;
                for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                    if ( i != xStartCell){
                        cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                        pR[cellNumber][pR_size[cellNumber]] = r;
                        pR_size[cellNumber] = pR_size[cellNumber] + 1;
                    }
                    for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                        cellNumber = getCellId(i, j, runNumPartitionsPerRelation);

                        pR[cellNumber][pR_size[cellNumber]] = r;
                        pR_size[cellNumber] = pR_size[cellNumber] + 1;
                    }
                }
            }        
        }


        void insert(Record &r, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation ){

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            
            xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
            yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
            firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

            // Determine cell for (rec.xEnd, rec.yEnd)
            auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

            if (r.xEnd + EPS >= 1) {
                xEndCell = runNumPartitionsPerRelation - 1;
            }
            else {
                xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
            }

            if (r.yEnd + EPS >= 1) {
                yEndCell = runNumPartitionsPerRelation - 1;
            }
            else {
                yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
            }
            lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

            int x = 0 , y = 0;

            // Put record in cells.
            if (firstCell == lastCell) {
                pRA[firstCell][pRA_size[firstCell]] = r;
                pRA_size[firstCell] = pRA_size[firstCell] + 1;
            }
            else {
                pRA[firstCell][pRA_size[firstCell]] = r;
                pRA_size[firstCell] = pRA_size[firstCell] + 1;
                
                int cellNumber;
                for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                    if ( i != xStartCell){
                        cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                        pRC[cellNumber][pRC_size[cellNumber]] = r;
                        pRC_size[cellNumber] = pRC_size[cellNumber] + 1;
                    }
                    for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                        cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                        if( i == xStartCell){
                            pRB[cellNumber][pRB_size[cellNumber]] = r;
                            pRB_size[cellNumber] = pRB_size[cellNumber] + 1 ;
                        }
                        else{
                            pRD[cellNumber][pRD_size[cellNumber]] = r;
                            pRD_size[cellNumber] = pRD_size[cellNumber] + 1 ;
                        }
                    }
                }
            }
        }


        void Partition(const Relation& R, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;
            Timer tim;
            double timepR = 0, timeDecomp = 0;

            for (size_t i = 0; i < R.numRecords; i++){
                auto &r = R[i];

                // Determine cell for (rec.xStart, rec.yStart)
                xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }

                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }

                lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                int x = 0, y = 0;

                // Put record in cells.
                if (firstCell == lastCell) {
                    pRA_size[firstCell]++;
                }
                else {
                    pRA_size[firstCell]++;
                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                            pRC_size[cellNumber]++;
                        }
                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                            if( i == xStartCell){
                                pRB_size[cellNumber]++;
                            }
                            else{
                                pRD_size[cellNumber]++;
                            }
                        }
                    }
                }
            }

            int counter = 0;
            for (int i = 0; i < runNumPartitions; i++){
                pRA[i].resize(pRA_size[i]);
                pRB[i].resize(pRB_size[i]);
                pRC[i].resize(pRC_size[i]);
                pRD[i].resize(pRD_size[i]);                        
            }
            
            for (int i = 0; i < runNumPartitions; i++){
                pRD_size[i] = 0;
                pRC_size[i] = 0;
                pRB_size[i] = 0;
                pRA_size[i] = 0;           
                
            }

            for (size_t i = 0; i < R.numRecords; i++){
                auto &r = R[i];

                xStartCell = myQuotient(r.xStart + EPS, partitionExtent);
                yStartCell = myQuotient(r.yStart + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }

                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }
                lastCell = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);

                int x = 0 , y = 0;

                // Put record in cells.
                if (firstCell == lastCell) {
                    pRA[firstCell][pRA_size[firstCell]] = r;
                    
                    pRA_size[firstCell] = pRA_size[firstCell] + 1;   
                }
                else {
                    pRA[firstCell][pRA_size[firstCell]] = r;
                                
                    pRA_size[firstCell] = pRA_size[firstCell] + 1;

                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);

                            pRC[cellNumber][pRC_size[cellNumber]] = r;
                            
                            pRC_size[cellNumber] = pRC_size[cellNumber] + 1;

                        }

                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                            if( i == xStartCell){

                                pRB[cellNumber][pRB_size[cellNumber]] = r;
                                
                                pRB_size[cellNumber] = pRB_size[cellNumber] + 1 ;

                            }
                            else{

                                pRD[cellNumber][pRD_size[cellNumber]] = r;
                                
                                pRD_size[cellNumber] = pRD_size[cellNumber] + 1 ;
                            }
                        }
                    }
                }
            }
        };


        void PartitionTwoDimensional(Relation& R, Relation *pRA, Relation *pRB, Relation *pRC, Relation *pRD, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size, size_t *pRD_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
            Partition(R, pRA, pRB, pRC, pRD, pRA_size, pRB_size, pRC_size, pRD_size, runNumPartitionsPerRelation);            
        };

         void PartitionUniform(const Relation& R, Relation *pR, size_t *pR_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            

            int vec_extend = 0;

            double xStartCell, yStartCell, xEndCell, yEndCell;
            int firstCell, lastCell;

            for (size_t i = 0; i < R.numRecords; i++){
                auto &r = R[i];
                // Determine cell for (rec.xStart, rec.yStart)
                xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
                yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }
                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }


                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation-1;
                }
                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }

                lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

                // Put record in cells.
                if (firstCell == lastCell) {
                    pR_size[firstCell]++;

                }
                else {
                    pR_size[firstCell]++;

                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                            pR_size[cellNumber]++;
                        }
                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);

                            pR_size[cellNumber]++;
                        }
                    }
                }
            }

            for (int i = 0; i < runNumPartitions; i++){
                pR[i].resize(pR_size[i]);
                pR_size[i] = 0;
            }

            for (size_t i = 0; i < R.numRecords; i++){
                
                auto &r = R[i];

                xStartCell = myQuotient(r.xStart + EPS , partitionExtent);
                yStartCell = myQuotient(r.yStart  + EPS, partitionExtent);
                firstCell = getCellId(xStartCell, yStartCell,runNumPartitionsPerRelation);

                // Determine cell for (rec.xEnd, rec.yEnd)
                auto xEnd = myRemainder(r.xEnd, partitionExtent, int(myQuotient(r.xEnd + EPS, partitionExtent)));
                auto yEnd = myRemainder(r.yEnd, partitionExtent, int(myQuotient(r.yEnd + EPS, partitionExtent)));

                if (r.xEnd + EPS >= 1) {
                    xEndCell = runNumPartitionsPerRelation - 1;
                }

                else {
                    xEndCell = myQuotient(r.xEnd + EPS, partitionExtent);
                }

                if (r.yEnd + EPS >= 1) {
                    yEndCell = runNumPartitionsPerRelation-1;
                }

                else {
                    yEndCell = myQuotient(r.yEnd + EPS, partitionExtent);
                }

                lastCell = getCellId(xEndCell, yEndCell,runNumPartitionsPerRelation);

                // Put record in cells.
                if (firstCell == lastCell) {
                    pR[firstCell][pR_size[firstCell]] = r ;
                    pR_size[firstCell]++;
                }
                else {
                    pR[firstCell][pR_size[firstCell]] = r ;
                    pR_size[firstCell]++;

                    int cellNumber;
                    for ( int i = xStartCell ; i <= xEndCell ; i++ ){
                        if ( i != xStartCell){
                            cellNumber = getCellId(i, yStartCell, runNumPartitionsPerRelation);
                            pR[cellNumber][pR_size[cellNumber]] = r ;
                            pR_size[cellNumber]++;
                            
                        }
                        for ( int j = yStartCell + 1 ; j <= yEndCell ; j ++ ){
                            cellNumber = getCellId(i, j, runNumPartitionsPerRelation);
                            pR[cellNumber][pR_size[cellNumber]] = r ;
                            pR_size[cellNumber]++;
                        }
                    }
                }        
            }
        }
            

        void PartitionTwoDimensional(const Relation& R, Relation *pR, size_t *pR_size, int runNumPartitionsPerRelation)
        {
            int runNumPartitions = runNumPartitionsPerRelation * runNumPartitionsPerRelation;
            PartitionUniform(R, pR, pR_size, runNumPartitionsPerRelation);
        };
    }

}

