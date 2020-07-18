#ifndef TWOLEVELPLUS_H
#define	TWOLEVELPLUS_H

namespace twoLevelPlus{ 
    
    int myRemainder(int numer, int denom) {
        return numer%denom;
    };
    
    double myRemainder2(double numer, double denom, int q) {
        double rem = double(numer - q*denom);

        return ((abs(rem) < EPS) ? 0: rem);
    };


    int myQuotient(double numer, double denom) {
        return int(numer/denom + EPS);
    };
    

    int findReferenceCell(double x, double y, double cellExtent, int numCellsPerDimension) {
        int xInt,yInt;

        xInt = (x + EPS)/cellExtent;
        yInt = (y + EPS)/cellExtent;

        return (yInt * numCellsPerDimension + xInt);
    };

    int getCellId(int x, int y, int numCellsPerDimension) {
        return (y * numCellsPerDimension + x);
    };

    int FindRelevantTiles(Record &s,vector<int> &cornerCells,vector<int> &insideCells,int runNumPartitionsPerRelation)
    {
        double xStartCell, yStartCell, xEndCell, yEndCell;

        int queryCase;
        int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
        Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
        xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
        yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
        auto xEnd = myRemainder2(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
        auto yEnd = myRemainder2(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));
        

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
    
    namespace window
    {
        namespace sort{

                    void SortXStart(vector<Decompose> *pR, size_t *pRA_size,  int runNumPartitions){
        
                        for (int i = 0; i < runNumPartitions; i++){
                            std::sort(pR[i].begin(), pR[i].begin() + pRA_size[i]);
                        }
                    };

                    void SortXEnd(vector<Decompose> *pR, size_t *pRA_size, size_t *pRB_size, size_t *pRC_size,  int runNumPartitions){
        
                        for (int i = 0; i < runNumPartitions; i++){
                            std::sort(pR[i].begin(), pR[i].begin() + pRA_size[i]);
                            std::sort(pR[i].begin() + pRB_size[i], pR[i].begin() + pRC_size[i]);
                        }
                    };

                    void SortYStart(vector<Decompose> *pR, size_t *pRA_size,  int runNumPartitions){
        
                        for (int i = 0; i < runNumPartitions; i++){
                            std::sort(pR[i].begin(), pR[i].begin() + pRA_size[i]);
                        }
                    };

                    void SortYEnd(vector<Decompose> *pR, size_t *pRA_size, size_t *pRB_size,  int runNumPartitions){
        
                        for (int i = 0; i < runNumPartitions; i++){
                            std::sort(pR[i].begin(), pR[i].begin() + pRA_size[i]);
                            std::sort(pR[i].begin() + pRA_size[i], pR[i].begin() + pRB_size[i]);
                        }
                    };

            };

            inline void InternalLoop_Range(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
            {                 
                auto pivot = firstFS;
                while ((pivot < lastFS))
                {
                    if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart || S.xStart > pivot->xEnd || S.xEnd < pivot->xStart)
                    {
                        pivot ++;
                        continue;
                    }
                    resultItemsIds.push_back(pivot->id);
                    pivot ++;
                }
            }

            inline void RangeCorners(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
            {
                auto r = pR.begin() + startR;
                auto lastR = pR.begin() + endR;

                InternalLoop_Range(r,S,lastR,resultItemsIds);
            };

            inline void InternalLoop_Range_B(vector<Decompose>::const_iterator firstFS, vector<Decompose>::const_iterator lastFS, Record &s, vector<double> &indexFS,vector<RecordId> &resultItemsIds)
            { 
                auto pivot = firstFS;
                while ((pivot < lastFS))
                {
                    if ( (s.yStart > indexFS[(*pivot).id*4+3]) || (s.yEnd < indexFS[(*pivot).id*4+2]))
                    {
                        pivot++;
                        continue;
                    }
                    resultItemsIds.push_back(pivot->id);
                    pivot++;
                }
            }

            inline void RangeBClass(vector<Decompose> &pR,vector<Coord> &indexR, Record &s, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
            {
                auto r = pR.begin() + startR ;
                auto lastR = pR.begin() + endR;

                InternalLoop_Range_B(r, lastR/*, count*/,s, indexR,resultItemsIds);

            };

            inline void InternalLoop_Range_C(vector<Decompose>::const_iterator firstFS, vector<Decompose>::const_iterator lastFS, Record &s, vector<double> &indexFS,vector<RecordId> &resultItemsIds)
            { 
                auto pivot = firstFS;
                while ((pivot < lastFS))
                {
                    if ((s.xStart > indexFS[(*pivot).id*4+1]) || (s.xEnd < indexFS[(*pivot).id*4]))
                    {
                        pivot++;
                        continue;
                    }
                    resultItemsIds.push_back(pivot->id);
                    pivot++;
                }
            }

            inline void RangeCClass(vector<Decompose> &pR,vector<Coord> &indexR, Record &s, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
            {
                auto r = pR.begin() + startR;
                auto lastR = pR.begin() + endR;

                InternalLoop_Range_C(r, lastR/*, count*/,s, indexR,resultItemsIds);
            };

            int binarySearch(vector<Decompose> &R, size_t startR, size_t endR,double key)
            {
                if(startR >= endR)
                    return 0;
                
                int low = startR;
                int high = endR;
                
                int ans = endR;
                bool flag = false;
                while (low <= high) {
                    
                    int mid = (high+low)/2;
                    double midVal = R[mid].value;
                    
                    if (midVal < key) {
                        low = mid + 1;
                    }
                    else if (midVal >= key) {
                        ans = mid;
                        high = mid - 1;
                    }
                }
                if(low > endR)
                    return 0;
                return endR-ans;
            }

            inline void RangeWithBinarySearch(vector<Decompose> &pR, size_t startR, size_t endR,double key,vector<RecordId> &resultItemsIds)
            {
                int index = binarySearch(pR, startR, endR,key);
                
                if (index > 0)
                {
                    for(int i = endR-index; i < endR; i++)
                    {
                        resultItemsIds.push_back(pR[i].id);
                    }
                }
            };  

            int binarySearch2(vector<Decompose> &R, size_t startR, size_t endR,double key)
            {
                if(startR >= endR)
                    return 0;
                
                int low = startR;
                int high = endR;
                
                int ans = endR;

                while (low <= high) {
                    
                    int mid = (high+low)/2;
                    double midVal = R[mid].value;
                    
                    if (midVal < key) {
                        low = mid + 1;
                    }
                    else if (midVal > key) {
                        ans = mid;
                        high = mid - 1;
                    }
                    else if (midVal == key) {
                        low = mid+1;
                    }
                }
                if(low > endR){
                    return endR;
                }
                
                return ans;
            }

            inline void RangeWithBinarySearch2(vector<Decompose> &pR, size_t startR, size_t endR,double key,vector<RecordId> &resultItemsIds)
            {
                int index = binarySearch2(pR, startR, endR,key);
                if (index >= 0)
                {
                    for(int i = 0; i < index;i++)
                    {
                        resultItemsIds.push_back(pR[i].id);
                    }
                }
            }; 

    }

    namespace disk
    {
        
    }
}

#endif

