#ifndef BATCHPROCESSING_H
#define	BATCHPROCESSING_H


namespace batchProcessing
{
    namespace window
    {
        namespace vec
        {

            inline void InternalLoop_Rolled_Range(RelationIterator firstFS,Record &S, RelationIterator lastFS, vector<RecordId> &resultItemsIds)
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

            inline void RangeCorners(Relation &R, Record &S, size_t endR, vector<RecordId> &resultItemsIds)
            {
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                InternalLoop_Rolled_Range(r,S,lastR,resultItemsIds);
            };

            inline void InternalLoop_Rolled_Range_B_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS, vector<RecordId> &resultItemsIds)
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
            };

            inline void RangeBClass(Relation &R, Record &S, size_t endR, vector<RecordId> &resultItemsIds)
            {
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                InternalLoop_Rolled_Range_B_Class(r,S,lastR,resultItemsIds);
            };

            inline void InternalLoop_Rolled_Range_C_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS, vector<RecordId> &resultItemsIds)
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
            };

            inline void RangeCClass(Relation &R, Record &S, size_t endR, vector<RecordId> &resultItemsIds)
            {
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                InternalLoop_Rolled_Range_C_Class(r,S,lastR,resultItemsIds);
            };
        };

        namespace xor_workload
        {

            inline double InternalLoop_Rolled_Range(RelationIterator firstFS,Record &S, RelationIterator lastFS)
            { 
                unsigned long long result = 0;
                
                auto pivot = firstFS;
                while ((pivot < lastFS))
                {
                    if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart || S.xStart > pivot->xEnd || S.xEnd < pivot->xStart)
                    {
                        pivot ++;
                        continue;
                    }
                    result += pivot->id ^ S.id;
                    pivot ++;
                }
                return result;
            }

            inline unsigned long long RangeCorners(Relation &R, Record &S, size_t endR)
            {
                unsigned long long result = 0;
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                result += InternalLoop_Rolled_Range(r,S,lastR);

                return result;
            };

            inline unsigned long long InternalLoop_Rolled_Range_B_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS)
                { 
                    unsigned long long result = 0;
                    
                    auto pivot = firstFS;
                    while ((pivot < lastFS))
                    {    
                    if (S.yStart > pivot->yEnd || S.yEnd < pivot->yStart )
                        {
                            pivot++;
                            continue;
                        }
                        result+= pivot->id ^ S.id;
                        pivot++;
                    }
                    return result;
                }

            inline unsigned long long RangeBClass(Relation &R, Record &S, size_t endR)
            {
                unsigned long long result = 0;
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                result += InternalLoop_Rolled_Range_B_Class(r,S,lastR);

                return result;
            };

            inline unsigned long long InternalLoop_Rolled_Range_C_Class(RelationIterator firstFS,Record &S, RelationIterator lastFS)
            { 
                unsigned long long result = 0;
                auto pivot = firstFS;
                while ((pivot < lastFS))
                {
                    if ( S.xStart > pivot->xEnd || S.xEnd < pivot->xStart )
                    {
                        pivot++;
                        continue;
                    }
                    result+= pivot->id ^ S.id;
                    pivot++;
                }
                return result;
            }

            inline unsigned long long RangeCClass(Relation &R, Record &S, size_t endR)
            {
                unsigned long long result = 0;
                auto r = R.begin();
                auto lastR = R.begin() + endR;

                result += InternalLoop_Rolled_Range_C_Class(r,S,lastR);

                return result;
            };
        }
    }
}



#endif