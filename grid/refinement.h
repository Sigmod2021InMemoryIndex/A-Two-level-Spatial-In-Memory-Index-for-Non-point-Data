#ifndef REFINEMENT_H
#define	REFINEMENT_H

namespace refinement{

    Coord distBetweenPoints(Coordinates point, Coordinates centerPoint)
    {

        auto diffx = centerPoint.x - point.x;
        auto diffy = centerPoint.y - point.y;

        auto dist = sqrt((diffx*diffx) + (diffy*diffy));
        return dist;
    }

    unsigned long long refine(const Record& s, const Geometry &polygon)
    {
        size_t size = polygon.geometry.size();
        for(int i = 0; i < size; i ++)
        {
            Coordinates point = polygon.geometry.at(i);
            if(point.x >= s.xStart && point.x <= s.xEnd && point.y >= s.yStart && point.y <= s.yEnd)
            {
                return 1;
                break;
            }
        }
        return 0;
    }

    unsigned long long refineDisc(Coordinates centerPoint, const Geometry &polygon, Coord epsilon)
    {
        for(int i = 0; i < polygon.geometry.size(); i ++)
        {
            Coordinates point = polygon.geometry.at(i);
            if(distBetweenPoints(point,centerPoint) <= epsilon)
            {
                return 1;
            }
        }
        return 0;
    }

    void normalizeGeometries(Coord minX, Coord minY, Coord maxExtent,vector<Geometry> &geo)
    {
        for(auto &item : geo)
        {
            for(auto &point : item.geometry)
            {
                point.x = (point.x - minX) / maxExtent;
                point.y = (point.y - minY) / maxExtent;
            }
        }
    }


}







#endif