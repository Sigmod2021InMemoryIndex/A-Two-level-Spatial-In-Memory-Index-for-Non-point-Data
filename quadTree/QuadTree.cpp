#include <iostream>
#include "QuadTree.h"

BoundingBox::BoundingBox(double xStart, double yStart, double xEnd, double yEnd){
    this->xStart = xStart;
    this->yStart = yStart;
    this->xEnd = xEnd;
    this->yEnd = yEnd;
}


BoundingBox::BoundingBox(){
}


BoundingBox::~BoundingBox(){
}


QuadTreeNode::QuadTreeNode(BoundingBox boundary, int capacity, int level, int QUAD_TREE_SIZE){
    this->northWest = NULL;
    this->northEast = NULL;
    this->southWest = NULL;
    this->southEast = NULL;
    this->capacity = capacity;
    this->boundary = boundary;
    this->level = level;
    this->QUAD_TREE_SIZE = QUAD_TREE_SIZE;
    this->counter = 0; 
}


QuadTreeNode::QuadTreeNode(BoundingBox boundary, int capacity, int level){
    this->northWest = NULL;
    this->northEast = NULL;
    this->southWest = NULL;
    this->southEast = NULL;
    this->capacity = capacity;
    this->boundary = boundary;
    this->level = level;
    this->counter = 0; 
}


QuadTreeNode::~QuadTreeNode(){
}


void QuadTreeNode::loadData(Relation &R){
    
    for (size_t i = 0; i < R.numRecords; i++){
        auto &r = R[i];
        this->insert(r);
    }
}


void QuadTreeNode::loadDataNoHeight(Relation &R){
    
    for (size_t i = 0; i < R.numRecords; i++){
        auto &r = R[i];
        this->insertNoHeight(r);
    }
}


void QuadTreeNode::insert(Record r){
    if (inBoundary(r)){
        return;
    }
    
    if (counter < this->capacity){ //insert
        this->data.push_back(r);
        counter ++;
    }
    else if (this->level < this->QUAD_TREE_SIZE){
        if ( this->northWest == NULL ) {
            this->subdivide();
            
        }

        for ( int i = 0 ; i < this->data.size() ; i ++ ){
            auto &rec =  this->data.at(i);
            this->northWest->insert (rec);
            this->northEast->insert (rec);
            this->southWest->insert (rec);
            this->southEast->insert (rec);
        }
        
        this->data.clear();

        this->northWest->insert (r);
        this->northEast->insert (r);
        this->southWest->insert (r);
        this->southEast->insert (r);

    }
    else{ // subdivide     
        this->data.push_back(r);
        counter ++;
    }
}


void QuadTreeNode::insertNoHeight(Record r){
    if (inBoundary(r)){
        return;
    }

    if (counter < this->capacity){ //insert
        this->data.push_back(r);
        counter ++;
    }
    else{
        if ( this->northWest == NULL ) {
            this->subdivideNoHeight();     
        }

        for ( int i = 0 ; i < this->data.size() ; i ++ ){
            auto &rec =  this->data.at(i);
            this->northWest->insertNoHeight (rec);
            this->northEast->insertNoHeight (rec);
            this->southWest->insertNoHeight (rec);
            this->southEast->insertNoHeight (rec);
        }
        
        this->data.clear();

        this->northWest->insertNoHeight (r);
        this->northEast->insertNoHeight (r);
        this->southWest->insertNoHeight (r);
        this->southEast->insertNoHeight (r);
    }
}


unsigned long long QuadTreeNode::calculateSize(){
    unsigned long long s = 0;

    if ( this->northWest == NULL ) {
        s = this->data.size()*(4*sizeof(Coord)+sizeof(RecordId));

        return s;
    }

    s += this->northWest->calculateSize();
    s += this->northEast->calculateSize();
    s += this->southWest->calculateSize();
    s += this->southEast->calculateSize();
    
    return s;
}


void QuadTreeNode::deleteTree(){

    if ( this->northWest == NULL ) {
        this->northWest->deleteTree();
        this->northEast->deleteTree();
        this->southWest->deleteTree();
        this->southEast->deleteTree();
    }

    delete this;
}


void QuadTreeNode::subdivideNoHeight(){
    BoundingBox box = this->boundary;
    
    // Spit the quadrant into four equal parts.
    double xMid = (box.xEnd + box.xStart) / 2.0;
    double yMid = (box.yEnd + box.yStart) / 2.0;

    // Create the north west bounding box.
    BoundingBox northWest = BoundingBox( box.xStart, yMid, xMid, box.yEnd );
    this->northWest = new QuadTreeNode ( northWest, this->capacity,this->level + 1);
    
    // Create the north east bounding box.
    BoundingBox northEast = BoundingBox ( xMid, yMid, box.xEnd, box.yEnd );
    this->northEast = new QuadTreeNode ( northEast, this->capacity, this->level + 1);
    
    // Create the south west bounding box.
    BoundingBox southWest = BoundingBox ( box.xStart, box.yStart, xMid, yMid );
    this->southWest = new QuadTreeNode ( southWest, this->capacity, this->level+1);
    
    // Create the south east bounding box.
    BoundingBox southEast = BoundingBox ( xMid, box.yStart, box.xEnd, yMid );
    this->southEast = new QuadTreeNode ( southEast, this->capacity, this->level+1);
}



void QuadTreeNode::subdivide(){
    BoundingBox box = this->boundary;
    
    // Spit the quadrant into four equal parts.
    double xMid = (box.xEnd + box.xStart) / 2.0;
    double yMid = (box.yEnd + box.yStart) / 2.0;

    // Create the north west bounding box.
    BoundingBox northWest = BoundingBox( box.xStart, yMid, xMid, box.yEnd );
    this->northWest = new QuadTreeNode ( northWest, this->capacity,this->level + 1, this->QUAD_TREE_SIZE);
    
    // Create the north east bounding box.
    BoundingBox northEast = BoundingBox ( xMid, yMid, box.xEnd, box.yEnd );
    this->northEast = new QuadTreeNode ( northEast, this->capacity, this->level + 1, this->QUAD_TREE_SIZE);
    
    // Create the south west bounding box.
    BoundingBox southWest = BoundingBox ( box.xStart, box.yStart, xMid, yMid );
    this->southWest = new QuadTreeNode ( southWest, this->capacity, this->level+1, this->QUAD_TREE_SIZE);
    
    // Create the south east bounding box.
    BoundingBox southEast = BoundingBox ( xMid, box.yStart, box.xEnd, yMid );
    this->southEast = new QuadTreeNode ( southEast, this->capacity, this->level+1, this->QUAD_TREE_SIZE);
}


bool QuadTreeNode::inBoundary(Record r){
    return (this->boundary.xStart > r.xEnd || this->boundary.xEnd < r.xStart || this->boundary.yStart > r.yEnd || this->boundary.yEnd < r.yStart );
}


void QuadTreeNode::rangeQuery(Record s, vector<RecordId> &resultItemsIds){
    unsigned long long result = 0;

    if (inBoundary(s)){
        return;
    }

    if ( this->northWest == NULL ) {
        if ( s.xEnd > this->boundary.xStart || s.xStart < this->boundary.xEnd || s.yEnd  > this->boundary.yStart || s.yStart  < this->boundary.yEnd ){

            for ( int i = 0 ; i < this->data.size(); i ++){
                auto &r = this->data.at(i);

                if ((r.xStart > s.xEnd) || (r.xEnd < s.xStart) || (r.yStart > s.yEnd) || (r.yEnd < s.yStart)){
                    continue;
                }
                else{
                    auto x = max(s.xStart, r.xStart);
                    auto y = min(s.yEnd, r.yEnd);
                    
                    if (this->boundary.xStart > x || x > this->boundary.xEnd || this->boundary.yStart > y || y > this->boundary.yEnd){
                        continue;
                    }
                    resultItemsIds.push_back(r.id);
                }

            }
        }
        else{
            for ( int i = 0 ; i < this->data.size(); i ++){
                auto &r = this->data.at(i);

                auto x = max(s.xStart, r.xStart);
                auto y = min(s.yEnd, r.yEnd);
                
                if (this->boundary.xStart > x || x > this->boundary.xEnd || this->boundary.yStart > y || y > this->boundary.yEnd){
                    continue;
                }
                resultItemsIds.push_back(r.id);
            }
        }
        return;
    }

    this->northWest->rangeQuery (s, resultItemsIds);
    this->northEast->rangeQuery (s, resultItemsIds);
    this->southWest->rangeQuery (s, resultItemsIds);
    this->southEast->rangeQuery (s, resultItemsIds);
}

Coord QuadTreeNode::MinDist(Coordinates &point, Record &rec){
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

};


Coord QuadTreeNode::MaxDist(Coordinates &point){
    Coord sum = 0.0;
    Coord diff1, diff2, diff;
    
    diff1 = abs(point.x - this->boundary.xStart);
    diff2 = abs(point.x - this->boundary.xEnd);
    
    diff = diff1 > diff2 ? diff1 : diff2; 
    sum = diff*diff;
              
    diff1 = 0.0;
    diff2 = 0.0;
    diff = 0.0;
       
    diff1 = abs(point.y - this->boundary.yStart);
    diff2 = abs(point.y - this->boundary.yEnd);
    
    diff = diff1 > diff2 ? diff1 : diff2; 
    
    sum += diff*diff; 
            
    return sum;
};


void QuadTreeNode::rangeQuery(Record s, vector<RecordId> &resultItemsIds, Coord epsilon){
    unsigned long long result = 0;

    if (inBoundary(s)){
        //return 0;
        return;
    }

    if ( this->northWest == NULL ) {

        Coord xCircle,yCircle;
    
        xCircle = (s.xStart + s.xEnd)/2;
        yCircle = (s.yStart + s.yEnd)/2;
        
        Coordinates point (xCircle,yCircle);

        if (MaxDist(point) <= epsilon*epsilon ){
            for ( int i = 0 ; i < this->data.size(); i ++){
                auto &r = this->data.at(i);

                auto x = max(s.xStart, r.xStart);
                auto y = min(s.yEnd, r.yEnd);
                
                if (this->boundary.xStart > x || x > this->boundary.xEnd || this->boundary.yStart > y || y > this->boundary.yEnd){
                    continue;
                }

                if ( MinDist(point,r) > epsilon*epsilon ){
                    continue;
                }

                resultItemsIds.push_back(r.id);    
            }
        }
        else{
            for ( int i = 0 ; i < this->data.size(); i ++){
                auto &r = this->data.at(i);

                if ((r.xStart > s.xEnd) || (r.xEnd < s.xStart) || (r.yStart > s.yEnd) || (r.yEnd < s.yStart)){
                    continue;
                }
                else{
                    auto x = max(s.xStart, r.xStart);
                    auto y = min(s.yEnd, r.yEnd);
                    
                    if (this->boundary.xStart > x || x > this->boundary.xEnd || this->boundary.yStart > y || y > this->boundary.yEnd){
                        continue;
                    }

                    if ( MinDist(point,r) > epsilon*epsilon ){
                        continue;
                    }

                    resultItemsIds.push_back(r.id);
                }

            }
        }
        return;
    }

    this->northWest->rangeQuery (s, resultItemsIds, epsilon);
    this->northEast->rangeQuery (s, resultItemsIds, epsilon);
    this->southWest->rangeQuery (s, resultItemsIds, epsilon);
    this->southEast->rangeQuery (s, resultItemsIds, epsilon);

};