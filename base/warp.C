#include "warp.h"
#include <iostream>

using namespace img;
using namespace std;

void img::ApplyWarp(const img::Point &cen, const double range,
                    const img::Warp &w, const img::ImgProc &in, img::ImgProc &out)
{


    out.clear(in.nx(), in.ny(), in.depth());
    for(int j = 0; j < out.ny();j++) {
#pragma omp parallel for
        for(int i = 0; i < out.nx(); i++) {

        
        Point P;
        P.x = (float)i/(float)out.nx();
        P.y = (float)j/(float)out.ny();
        P.x *= range;
        P.y *= range;

        P.x -= cen.x;
        P.y -= cen.y;
        Point PP = w(P);
        PP.x += cen.x;
        PP.y += cen.y;
        PP.x /= range;
        PP.y /= range;
        PP.x *= in.nx();
        PP.y *= in.ny();
        vector<float> v(3);
 
        in.value_in(PP.x, PP.y,v);
        out.set_value(i,j,v);
        }
    }
}

Translate::Translate(const Point& t): trans(t) {}
Point Translate::operator()(const Point& P) const {
    Point result;
    result.x = P.x - trans.x;
    result.y = P.y - trans.y;
    return result;
}
