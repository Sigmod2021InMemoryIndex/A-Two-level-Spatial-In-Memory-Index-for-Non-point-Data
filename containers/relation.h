#ifndef _RELATION_H_
#define _RELATION_H_

#include "../def.h"

class Coordinates {
    public:
	Coord x,y;

	Coordinates(double x, double y);
        Coordinates();
	void print();
};


class Geometry {
    public:
	vector<Coordinates> geometry;
	Geometry(vector<Coordinates> geometry);
        Geometry();
	void print();
};


class Decompose {
    public: 
        RecordId id;
        Coord value;
        
        Decompose (RecordId id, Coord value);
        Decompose ();
        bool operator < (const Decompose& rhs) const;
        bool operator >= (const Decompose& rhs) const;
        bool operator < (const double rh) const;
};       

class Record {
    public:
        RecordId id;

        // MBR
        Coord xStart, yStart; // bottom-left corner
        Coord xEnd, yEnd;     // top-right corner

        Record();
        Record(RecordId id, Coord xStart, Coord yStart, Coord xEnd, Coord yEnd);
        bool operator < (const Record& rhs) const;
        bool operator >= (const Record& rhs) const;
        void print() const;
        void print(char c) const;
        ~Record();	
};



class Relation : public vector<Record>
{
    public:
	size_t numRecords;
        Coord minX, maxX, minY, maxY;
        double avgXExtent, avgYExtent;

	Relation();
        void loadDisk(Coord epsilon);
	void load(const char *filename);
        void loadWithGeometries(const char *filename,vector<Geometry> &geo);
	void sortByXStart();
        void sortByYStart();
        void normalize(Coord minX, Coord maxX, Coord minY, Coord maxY, Coord maxExtent);
        void computeAvgExtents1d();
        ~Relation();
};


typedef Relation::const_iterator RelationIterator;

#endif //_RELATION_H_
