#ifndef ONELEVEL_H
#define	ONELEVEL_H

namespace oneLevel{

    // int myQuotient(int numer, int denom) {
    //     return numer/denom;
    // };

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

    unsigned long long checkDuplicates(const Record& s,const Record& polygon,int cellNumber,int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
    {
        unsigned long long result = 0;
        Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
        double x,y;

        x = max(polygon.xStart, s.xStart);
        y = min(polygon.yEnd, s.yEnd);
        
        
        auto pid_ref = findReferenceCell(x, y, partitionExtent, runNumPartitionsPerRelation);

        if (pid_ref == cellNumber)
        {
            resultItemsIds.push_back(polygon.id);
        }
        return result;
    }

    namespace window
    {
        inline void RangeInside(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            for (int pid = 0; pid < size; pid++)
            { 
                checkDuplicates(s, pR[item][pid], item,runNumPartitionsPerRelation,resultItemsIds);
            }
        }

        inline void RangeCorner(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            for (int pid = 0; pid < size; pid++)
            {
                if (s.yStart > pR[item][pid].yEnd || s.yEnd < pR[item][pid].yStart || s.xStart > pR[item][pid].xEnd || s.xEnd < pR[item][pid].xStart)
                {
                    continue;
                }
                else
                {
                    checkDuplicates(s,pR[item][pid],item,runNumPartitionsPerRelation,resultItemsIds);                               
                }
            }
        }

        inline void RangeBottomBorder(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            for (int pid = 0; pid < size; pid++)
            {
                if (s.yStart > pR[item][pid].yEnd || s.yEnd < pR[item][pid].yStart)
                {
                    continue;
                }
                else
                {
                    checkDuplicates(s,pR[item][pid],item,runNumPartitionsPerRelation,resultItemsIds);                                 
                }
                
            }
        }

        inline void RangeLeftBorder(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            for (int pid = 0; pid < size; pid++)
            {
                if (s.xStart > pR[item][pid].xEnd || s.xEnd < pR[item][pid].xStart)
                {
                    continue;
                }
                else
                {
                    checkDuplicates(s,pR[item][pid],item,runNumPartitionsPerRelation,resultItemsIds);                                 
                }
                
            }
        }


        inline void RangeRightBorder(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            
            for (int pid = 0; pid < size; pid++)
            {
                if (s.xStart > pR[item][pid].xEnd || s.xEnd < pR[item][pid].xStart)
                {
                    continue;
                }
                else
                {
                    checkDuplicates(s,pR[item][pid],item,runNumPartitionsPerRelation,resultItemsIds);                                 
                }
                
            }
        }

        
        inline void RangeTopBorder(Relation *pR,int item, Record &s, int runNumPartitionsPerRelation,vector<RecordId> &resultItemsIds)
        {
            int size = pR[item].size();
            
            for (int pid = 0; pid < size; pid++)
            {
                if (s.yStart > pR[item][pid].yEnd || s.yEnd < pR[item][pid].yStart )
                {
                    continue;
                }
                else
                {
                    checkDuplicates(s,pR[item][pid],item,runNumPartitionsPerRelation,resultItemsIds);   
                }
                
            }
        }
    }

    namespace disk
    {
        Coord MaxDist(Coordinates &point, Record &rec, int dimensionality){
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

        void FindRelevantTiles(Record &s, Relation &cellsCoord, int runNumPartitionsPerRelation)
        {
            Coord xStartCell, yStartCell, xEndCell, yEndCell;
            int queryCase;
            int runNumPartitions = runNumPartitionsPerRelation*runNumPartitionsPerRelation;
            Coord partitionExtent = 1.0/runNumPartitionsPerRelation;
            xStartCell = myQuotient(s.xStart + EPS, partitionExtent);
            yStartCell = myQuotient(s.yStart + EPS, partitionExtent);
            auto xEnd = myRemainder2(s.xEnd , partitionExtent, int(myQuotient(s.xEnd + EPS, partitionExtent)));
            auto yEnd = myRemainder2(s.yEnd , partitionExtent, int(myQuotient(s.yEnd + EPS, partitionExtent)));
            
            

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
                yEndCell = myQuotient(s.yEnd  + EPS, partitionExtent);
            }

            int bottomLeft = getCellId(xStartCell, yStartCell, runNumPartitionsPerRelation);
            int bottomRight = getCellId(xEndCell, yStartCell, runNumPartitionsPerRelation);
            int topLeft = getCellId(xStartCell, yEndCell, runNumPartitionsPerRelation);
            int topRight = getCellId(xEndCell, yEndCell, runNumPartitionsPerRelation);
            
            int rightLimit = bottomRight;
            int start = bottomLeft;
            int ii = start; 
            while( ii<=topRight)
            {
                if(ii <= rightLimit)
                {

                    Coord xStart, yStart, xEnd, yEnd;
                    
                    xStart = myRemainder2(ii , runNumPartitionsPerRelation, int(myQuotient(ii ,runNumPartitionsPerRelation))) * partitionExtent;
                    yStart = myQuotient(ii, runNumPartitionsPerRelation) * partitionExtent ;

                    xEnd = xStart + partitionExtent;
                    yEnd = yStart + partitionExtent;
                    
                    cellsCoord.emplace_back(ii, xStart,yStart,xEnd, yEnd);
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
};

#endif