#ifndef TWOLEVEL_H
#define	TWOLEVEL_H

namespace twoLevel{ 
    
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
    
    namespace window
    {
    
        inline void InternalLoop_Range_Corners(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
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

        inline void Range_Corners(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_Corners(r,S,lastR,resultItemsIds);
        };

        inline void InternalLoop_Range_B_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
        {                 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {    
                if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart )
                {
                    pivot++;
                    continue;
                }
                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void Range_B_Class(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_B_Class(r,S,lastR,resultItemsIds);
        };


        inline void InternalLoop_Range_C_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
        { 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {
                if ( S.xStart > pivot->xEnd || S.xEnd < pivot->xStart )
                {
                    pivot++;
                    continue;
                }
                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void Range_C_Class(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_C_Class(r,S,lastR,resultItemsIds);

        };

    
    }

    namespace disk
    {
        void load(Relation &S, Coord epsilon){
            Coord x,y;
            
            for ( int i = 0 ; i < S.numRecords ; i ++){
                x = (S[i].xStart + S[i].xEnd)/2;
                y = (S[i].yStart + S[i].yEnd)/2;
                
                S[i].xStart = x - epsilon; 
                S[i].xEnd = x + epsilon;
                S[i].yStart = y - epsilon;
                S[i].yEnd = y + epsilon;
            }

        }

        Coord distPointToMbr(double x,double y, Coordinates centerPoint)
        {
            auto diffx = centerPoint.x - x;
            auto diffy = centerPoint.y - y;

            auto dist = sqrt((diffx*diffx) + (diffy*diffy));
            return dist;
        }

        void findCoordinates(int id , int &xStart, int &yStart, int runNumPartitionsPerRelation){
            xStart = id%runNumPartitionsPerRelation;
            yStart = id/runNumPartitionsPerRelation;
        }

        Coord MaxDist(Coordinates &point, Record &rec){
            Coord sum = 0.0;
            Coord diff1, diff2, diff;
            
            diff1 = abs(point.x - rec.xStart);
            diff2 = abs(point.x - rec.xEnd);
            
            diff = diff1 > diff2 ? diff1 : diff2; 
            sum = diff*diff;
                    
            
            diff1 = 0.0;
            diff2 = 0.0;
            diff = 0.0;
            
            
            diff1 = abs(point.y - rec.yStart);
            diff2 = abs(point.y - rec.yEnd);
            
            diff = diff1 > diff2 ? diff1 : diff2; 
            
            sum += diff*diff; 
                    
            return sum;
        }


        Coord MinDist(Coordinates &point, Record &rec){
            Coord sum = 0.0;
            Coord diff;

            diff = 0.0;
            if ( point.x < rec.xStart) {
                diff = rec.xStart - point.x;
            }
            else if ( point.x > rec.xEnd ){
                diff = point.x-rec.xEnd;
            }
            sum = diff*diff;
            diff = 0.0;

            if ( point.y < rec.yStart ) {
                diff = rec.yStart - point.y;
            }
            else if ( point.y > rec.yEnd ){
                diff = point.y - rec.yEnd;
            }
            sum += diff*diff;
                    
            return sum;

        }

        Coord MinDist(Coordinates &point, RelationIterator &rec){
            Coord sum = 0.0;
            Coord diff;

            diff = 0.0;
            if ( point.x < rec->xStart) {
                diff = rec->xStart - point.x;
            }
            else if ( point.x > rec->xEnd ){
                diff = point.x-rec->xEnd;
            }
            sum = diff*diff;
            diff = 0.0;

            if ( point.y < rec->yStart ) {
                diff = rec->yStart - point.y;
            }
            else if ( point.y > rec->yEnd ){
                diff = point.y - rec->yEnd;
            }
            sum += diff*diff;
                    
            return sum;

        }

        inline void InternalLoop_Range_Disk(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds, Coord epsilon, Coordinates point)
        {                 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {
                if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart || S.xStart > pivot->xEnd || S.xEnd < pivot->xStart)
                {
                    pivot ++;
                    continue;
                }

                if (MinDist(point,pivot) > epsilon*epsilon ){
                    pivot ++;
                    continue;
                }

                resultItemsIds.push_back(pivot->id);
                pivot ++;
            }
        }

        inline void RangeCornersDisk(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds,Coord epsilon, Coordinates point)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_Disk(r,S,lastR,resultItemsIds, epsilon, point);
        };

        inline void InternalLoop_Range_B_Disk(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds, Coord epsilon, Coordinates point)
        {                 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {    
                if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart )
                {
                    pivot++;
                    continue;
                }

                if (MinDist(point,pivot) > epsilon*epsilon ){
                    pivot++;
                    continue;
                }

                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void RangeBClassDisk(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds, Coord epsilon, Coordinates point)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_B_Disk(r,S,lastR,resultItemsIds, epsilon, point);
        };

        inline void InternalLoop_Range_C_Disk(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds, Coord epsilon, Coordinates point)
        { 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {
                if ( S.xStart > pivot->xEnd || S.xEnd < pivot->xStart )
                {
                    pivot++;
                    continue;
                }

                if (MinDist(point,pivot) > epsilon*epsilon ){
                    pivot++;
                    continue;
                }
                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void RangeCClassDisk(Relation &pR, Record &S, size_t startR, size_t endR,vector<RecordId> &resultItemsIds, Coord epsilon, Coordinates point)
        {
            auto r = pR.begin() + startR;
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_C_Disk(r,S,lastR,resultItemsIds, epsilon, point);
        };

        int FindRelevantTiles(Record &s,vector<int> &cornerCells,vector<int> &insideCells,int runNumPartitionsPerRelation)
        {
            double xStartCell, yStartCell, xEndCell, yEndCell;

            // cout << s.xStart << " " << s.xEnd << " " << s.yStart << " " << s.yEnd << endl;
            // cout << runNumPartitionsPerRelation << endl;

            int queryCase;
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
            yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
            auto xEnd = myRemainder2(s.xEnd, partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder2(s.yEnd, partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));
            

            if (s.xStart + EPS <= 0) {
                xStartCell = 0;
            }
            else {
                xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
            }
            
            if (s.yStart + EPS <= 0) {
                yStartCell = 0;
            }
            else {
                yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
            }
     
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
        
    }

    namespace update{

        inline void InternalLoop_Range_Corners(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
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

        inline void Range_Corners(Relation &pR, Record &S, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin();
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_Corners(r,S,lastR,resultItemsIds);
        };

        inline void InternalLoop_Range_B_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
        { 
            
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {    
                if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart )
                {
                    pivot++;
                    continue;
                }
                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void Range_B_Class(Relation &pR, Record &S, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin();
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_B_Class(r,S,lastR,resultItemsIds);
        };

        inline void InternalLoop_Range_C_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS,vector<RecordId> &resultItemsIds)
        { 
            auto pivot = firstFS;
            while ((pivot < lastFS))
            {
                if ( S.xStart > pivot->xEnd || S.xEnd < pivot->xStart )
                {
                    pivot++;
                    continue;
                }
                resultItemsIds.push_back(pivot->id);
                pivot++;
            }
        }

        inline void Range_C_Class(Relation &pR, Record &S, size_t endR,vector<RecordId> &resultItemsIds)
        {
            auto r = pR.begin();
            auto lastR = pR.begin() + endR;

            InternalLoop_Range_C_Class(r,S,lastR,resultItemsIds);

        };
    }
}

#endif

