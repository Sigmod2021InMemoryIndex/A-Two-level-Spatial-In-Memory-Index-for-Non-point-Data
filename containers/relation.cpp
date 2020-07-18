#include "relation.h"

bool CompareByYStart(const Record& lhs, const Record& rhs)
{
    return (lhs.yStart < rhs.yStart);
}


Coordinates::Coordinates( double x, double y)
{
    this->x = x;
    this->y = y;
}


Coordinates::Coordinates()
{

}


void Coordinates::print () {
    cout << this->x << " " << this->y << " \n";
}


Geometry::Geometry(vector<Coordinates> geometry){
    this->geometry = geometry;
}
        

Geometry::Geometry(){

}


Decompose::Decompose (RecordId id, Coord value){
    this->id = id;
    this->value = value;
}


Decompose::Decompose (){
}


bool Decompose::operator < (const Decompose& rhs) const
{
    return this->value < rhs.value;
}


bool Decompose::operator < (const double rh) const
{
    return this->value > rh;
}


bool Decompose::operator >= (const Decompose& rhs) const
{
    return !((*this) < rhs);
}


Record::Record()
{
}


Record::Record(RecordId id, Coord xStart, Coord yStart, Coord xEnd, Coord yEnd)
{
    this->id = id;

    // MBR
    this->xStart = xStart;
    this->xEnd   = xEnd;
    this->yStart = yStart;
    this->yEnd   = yEnd;
}


bool Record::operator < (const Record& rhs) const
{
    return this->xStart < rhs.xStart;
}

bool Record::operator >= (const Record& rhs) const
{
    return !((*this) < rhs);
}

   
Record::~Record()
{
}


Relation::Relation()
{
    this->minX = std::numeric_limits<Coord>::max();
    this->maxX = std::numeric_limits<Coord>::min();
    this->minY = std::numeric_limits<Coord>::max();
    this->maxY = std::numeric_limits<Coord>::min();
    this->avgXExtent = 0;
    this->avgYExtent = 0;
}


void Relation::load(const char *filename)
{
    RecordId id = 0;
    
    ifstream file( filename );
    if(!file)
    {
        cerr << "Cannot open the File : " << filename << endl;
        exit(1);
    }
    
    string line;
    while (getline( file, line ) )
    {
        
        istringstream is( line );
        string s;

        bool first = true;
        Coord minXmbr, minYmbr, maxXmbr, maxYmbr;
        minXmbr = std::numeric_limits<Coord>::max();
        maxXmbr = -std::numeric_limits<Coord>::max();
        minYmbr = std::numeric_limits<Coord>::max();
        maxYmbr = -std::numeric_limits<Coord>::max();
      
        //count vextor size
        int counter = 0 ;
        for ( int i = 0 ; i < line.size()  ; i ++){
        
            if ( line[i] == ','){
                counter++;
            }
        }
        counter++;

        while(getline( is, s, ',' ))
        {      
            Coord x,y;
            istringstream ss(s);
               
            ss >> x >> y;

            Coordinates c = Coordinates(x,y);
            
            this->minX = std::min(this->minX, x);
            this->maxX = std::max(this->maxX, x);
            this->minY = std::min(this->minY, y);
            this->maxY = std::max(this->maxY, y);

            minXmbr = std::min(minXmbr, x);
            maxXmbr = std::max(maxXmbr, x);
            minYmbr = std::min(minYmbr, y);
            maxYmbr = std::max(maxYmbr, y);
            
        }
        
        this->emplace_back(id, minXmbr, minYmbr, maxXmbr, maxYmbr);      
        id++;
    }

    file.close();

    this->numRecords = this->size();
}


void Relation::loadWithGeometries(const char *filename, vector<Geometry> &geo)
{
    RecordId id = 0;
    
    ifstream file( filename );
    if(!file)
    {
        cerr << "Cannot open the File : " << filename << endl;
        exit(1);
    }
    
    string line;
    while (getline( file, line ) )
    {       
        istringstream is( line );
        string s;

        bool first = true;
        Coord minXmbr, minYmbr, maxXmbr, maxYmbr;
        minXmbr = std::numeric_limits<Coord>::max();
        maxXmbr = -std::numeric_limits<Coord>::max();
        minYmbr = std::numeric_limits<Coord>::max();
        maxYmbr = -std::numeric_limits<Coord>::max();
      
        //count vextor size
        int counter = 0 ;
        for ( int i = 0 ; i < line.size()  ; i ++){
        
            if ( line[i] == ','){
                counter++;
            }
        }
        counter++;
        
        vector<Coordinates> vec ;      
        vec.reserve(counter);
        
        while(getline( is, s, ',' ))
        {            
            Coord x,y;
            istringstream ss(s);
               
            ss >> x >> y;

            Coordinates c = Coordinates(x,y);
            vec.push_back(c);
            
            this->minX = std::min(this->minX, x);
            this->maxX = std::max(this->maxX, x);
            this->minY = std::min(this->minY, y);
            this->maxY = std::max(this->maxY, y);

            minXmbr = std::min(minXmbr, x);
            maxXmbr = std::max(maxXmbr, x);
            minYmbr = std::min(minYmbr, y);
            maxYmbr = std::max(maxYmbr, y);
            
        }
        geo.push_back(vec);
        
        this->emplace_back(id, minXmbr, minYmbr, maxXmbr, maxYmbr);
        id++;
    }
    
    file.close();

    this->numRecords = this->size();
}


void Relation::loadDisk(Coord epsilon)
{
    Coord x,y;
    
    for ( int i = 0 ; i < this->numRecords ; i ++){
        x = (this->at(i).xStart +this->at(i).xEnd)/2;
        y = (this->at(i).yStart + this->at(i).yEnd)/2;
        
        this->at(i).xStart = x - epsilon; 
        this->at(i).xEnd = x + epsilon;
        this->at(i).yStart = y - epsilon;
        this->at(i).yEnd = y + epsilon;
    }
}


void Relation::normalize(Coord minX, Coord maxX, Coord minY, Coord maxY, Coord maxExtent) {
    for (Record& rec : (*this))
    { 
        rec.xStart = Coord(rec.xStart - minX) / maxExtent;
        rec.xEnd   = Coord(rec.xEnd   - minX) / maxExtent;
        rec.yStart = Coord(rec.yStart - minY) / maxExtent;
        rec.yEnd   = Coord(rec.yEnd   - minY) / maxExtent;
    }
}


void Relation::computeAvgExtents1d() {
    double sumX = 0, sumY = 0;

    for (Record& rec : (*this)) {
        sumX += rec.xEnd-rec.xStart;
        sumY += rec.yEnd-rec.yStart;
    }
    this->avgXExtent = (double)sumX/this->numRecords;
    this->avgYExtent = (double)sumY/this->numRecords;
}


Relation::~Relation()
{
}
