#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>
#include "./containers/relation.h"

using namespace std;

#define MAX_NODE_CAPACITY   16

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


void usage(){
    cerr << "NAME" << endl;
    cerr << "       ./rtree - range query using R-Tree" << endl << endl;
    cerr << "USAGE" << endl;
    cerr << "       ./rtree [OPTION]... [FILE1] [FILE2]" << endl << endl;
    cerr << "DESCRIPTION" << endl;
    cerr << "       -w" << endl;
    cerr << "              window query" << endl;
    cerr << "       -d" << endl;
    cerr << "              disk query. Option -e should be used when using -d " << endl;
    cerr << "       -e" << endl;
    cerr << "              radius of the disk query" << endl;
    cerr << "       Other arguments" << endl << endl;
    cerr << "       -i" << endl;
    cerr << "              number of iterations" << endl;
    cerr << "       -m" << endl;
    cerr << "              display this help message and exit" << endl << endl;
    cerr << "EXAMPLES" << endl;
    cerr << "       Window range query using R-Tree." << endl;
    cerr << "              ./rtree -w TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "       Disk range query using R-Tree with radius of 0.1." << endl;
    cerr << "              ./rtree -d -e 0.1 TIGER_ROADS_mbr.inp TIGER_c0.1%_n10000.qry" << endl;
    cerr << "\n" << endl;
    exit(1);
}


Coord MinDist(Coord pointX, Coord pointY, Coord xStartTile, Coord yStartTile, Coord xEndTile, Coord yEndTile, int dimensionality ){
    Coord sum = 0.0;
    Coord diff;

    diff = 0.0;
    if ( pointX < xStartTile) {
        diff = xStartTile - pointX;
    }
    else if ( pointX > xEndTile ){
        diff = pointX-xEndTile;
    }

    sum = diff*diff;
    diff = 0.0;

    if ( pointY < yStartTile ) {
        diff = yStartTile - pointY;
    }
    else if ( pointY > yEndTile ){
        diff = pointY - yEndTile;
    }
    
    sum += diff*diff;
               
    return sum;

}

int main(int argc, char **argv)
{
    char c;
    Relation R, S;
    double timeQuery = 0, indexBuild = 0;
    Timer tim;
    int nodeCapacity = -1, queryMethod=-1;
    Coord epsilon = -1.0;
    int NUM_ITERATIONS = 1;

    while ((c = getopt(argc, argv, "wdue:i:m")) != -1)
    {
        switch (c)
        {
            case 'w':
                queryMethod = WINDOW_QUERY;
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'd':
                queryMethod = DISK_QUERY;
                break;
            case 'i':
                NUM_ITERATIONS = atoi(optarg);
                break;
            case 'm':
                usage();
                break;
            default:
                cerr << "Wrong arguments! ";
                usage();
                break;
        }
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

    typedef bg::model::point<double, 2, bg::cs::cartesian> point;
    typedef bg::model::box<point> box;
    typedef std::pair<box, unsigned> value;
    
    if (queryMethod == DISK_QUERY){
        if (epsilon < 0){
            cout<<"Radius value is missing"<<endl;
            return 0;
        }
    }

    //Load inputs
    #pragma omp parallel sections
    {
         #pragma omp section
         {
            R.load(argv[optind]);
         }
         #pragma omp section
         {
            S.load(argv[optind+1]);
         }
    }

    cout << "Method\t" << "Query ID\t" << "Results\t" << "Indexing Time\t" << "Filtering Time" << endl;

    
    switch(queryMethod)
    {
        case WINDOW_QUERY:
        {
            Coord minX = min(R.minX, S.minX);
            Coord maxX = max(R.maxX, S.maxX);
            Coord minY = min(R.minY, S.minY);
            Coord maxY = max(R.maxY, S.maxY);
            Coord diffX = maxX - minX;
            Coord diffY = maxY - minY;
            Coord maxExtend = (diffX<diffY)?diffY:diffX;

            R.normalize(minX, maxX, minY, maxY, maxExtend);
            S.normalize(minX, maxX, minY, maxY, maxExtend);
            
            std::vector<box> contourCenters; // has some value
            std::vector<value> cloud;
        
            for ( unsigned i = 0 ; i < R.size() ; ++i )
            {
                box b(point(R[i].xStart, R[i].yStart), point(R[i].xEnd, R[i].yEnd));
                contourCenters.push_back(b);
            }

            tim.start();
            size_t id_gen = 0;
            std::transform(
                    contourCenters.begin(), contourCenters.end(),
                    back_inserter(cloud), 
                    [&](box const& p) { return std::make_pair(p, id_gen++); }
            );
            
            bgi::rtree<value, bgi::rstar<MAX_NODE_CAPACITY> > rtree(cloud.begin(),cloud.end());

            indexBuild = tim.stop(); 

            vector<box> queries;
            for ( unsigned i = 0 ; i < S.size() ; ++i )
            {
                box query_box(point(S[i].xStart, S[i].yStart), point(S[i].xEnd, S[i].yEnd));
                queries.push_back(query_box);
            }
            vector<double> indexTime[NUM_ITERATIONS];
            vector<double> queryTime[NUM_ITERATIONS];

            
            for ( int i = 0 ; i < queries.size() ; i ++){
                for ( int j = 0 ; j < NUM_ITERATIONS ; j++ ){
                    int results = 0;

                    tim.start();
                    
                    std::vector<value> result_s;
                    rtree.query(bgi::intersects(queries.at(i)), std::back_inserter(result_s));
                    
                    BOOST_FOREACH(value const& v, result_s)
                        results++;

                    timeQuery = tim.stop();
                    
                    std::cout<<"Rtree-Window-Query\t"<<i<<"\t"<<results<<"\t"<<indexBuild<<"\t"<< timeQuery<<endl;
                }
            }
            break;
        }

        case DISK_QUERY:
        {
            Coord minX = R.minX;
            Coord maxX = R.maxX;
            Coord minY = R.minY;
            Coord maxY = R.maxY;
            Coord diffX = maxX - minX;
            Coord diffY = maxY - minY;
            Coord maxExtend = (diffX<diffY)?diffY:diffX;

            R.normalize(minX, maxX, minY, maxY, maxExtend);
            S.normalize(minX, maxX, minY, maxY, maxExtend);

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
            
            std::vector<box> contourCenters; // has some value
            std::vector<value> cloud;
        
            for ( unsigned i = 0 ; i < R.size() ; ++i )
            {      
                box b(point(R[i].xStart, R[i].yStart), point(R[i].xEnd, R[i].yEnd));
                contourCenters.push_back(b);
            }

            tim.start();
            size_t id_gen = 0;
            std::transform(
                    contourCenters.begin(), contourCenters.end(),
                    back_inserter(cloud), 
                    [&](box const& p) { return std::make_pair(p, id_gen++); }
            );
            
            bgi::rtree<value, bgi::rstar<MAX_NODE_CAPACITY> > rtree(cloud.begin(),cloud.end());

            indexBuild = tim.stop(); 

            vector<box> queries;
            for ( unsigned i = 0 ; i < S.size() ; ++i )
            {
                
                box query_box(point(S[i].xStart, S[i].yStart), point(S[i].xEnd, S[i].yEnd));
                queries.push_back(query_box);
            }
            
            Coord xPoint;
            Coord yPoint;

            for ( int i = 0 ; i < queries.size() ; i ++){
                for(int o = 0; o < NUM_ITERATIONS; o++) {
                    tim.start();
                    auto &q = queries.at(i);
                    
                    
                    std::vector<value> returned_values;
                    point p = point((queries.at(i).min_corner().get<0>() + queries.at(i).max_corner().get<0>())/2,(queries.at(i).min_corner().get<1>() + queries.at(i).max_corner().get<1>())/2);
                    box b1(point(p.get<0>() - epsilon, p.get<1>() - epsilon), point(p.get<0>() + epsilon, p.get<1>() + epsilon));

                    xPoint = p.get<0>();  
                    yPoint = p.get<1>();
                    
                    rtree.query(bgi::intersects(b1),std::back_inserter(returned_values));

                    unsigned long long res = 0;
                    for (size_t i = 0; i < returned_values.size(); i++) {
                        double sum = 0.0;
                        double diff;
                        
                        if (MinDist(xPoint, yPoint, returned_values[i].first.min_corner().get<0>(), returned_values[i].first.min_corner().get<1>(), returned_values[i].first.max_corner().get<0>(), returned_values[i].first.max_corner().get<1>(), 2) <= epsilon*epsilon){
                            res++;
                        }
                        
                    }
                    timeQuery = tim.stop();

                    std::cout<<"Rtree-Disk-Query\t"<<i<<"\t"<<res<<"\t"<<indexBuild<<"\t"<< timeQuery<<endl;
                }
            }
            break;
        }
    }

    
    return 0;
}

