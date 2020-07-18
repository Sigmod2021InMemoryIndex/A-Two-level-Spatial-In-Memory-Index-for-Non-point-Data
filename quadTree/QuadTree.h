#ifndef _QuadTree_h_
#define _QuadTree_h_

#include "../containers/relation.h"

class BoundingBox{
    public :
        double xStart;
        double yStart;
        double xEnd;
        double yEnd;
        
        BoundingBox(double xStart, double yStart, double xEnd, double yEnd);
        BoundingBox();
        ~BoundingBox();
};

class QuadTreeNode{
    public :
        QuadTreeNode *northWest;
        QuadTreeNode *northEast;
        QuadTreeNode *southWest;
        QuadTreeNode *southEast;
        Relation data;

        int capacity;
        int id;
        int counter;
        int level;
        int QUAD_TREE_SIZE;
        BoundingBox boundary;

        QuadTreeNode(BoundingBox boundary, int capacity, int level, int QUAD_TREE_SIZE);
        QuadTreeNode(BoundingBox boundary, int capacity, int level);
        ~QuadTreeNode();
        void loadData(Relation &R);
        void loadDataNoHeight(Relation &R);
        void insert(Record r);
        void insertNoHeight(Record r);
        bool inBoundary(Record r); 
        void subdivide();
        void subdivideNoHeight();
        void rangeQuery(Record s, vector<RecordId> &resultItemsIds);
        void rangeQuery(Record s, vector<RecordId> &resultItemsIds, Coord epsilon);
        Coord MinDist(Coordinates &point, Record &rec);
        Coord MaxDist(Coordinates &point);
        void deleteTree();
        unsigned long long calculateSize();
};

#endif