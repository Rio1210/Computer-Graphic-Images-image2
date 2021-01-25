#include "imgproc.h"
using namespace std;
namespace img
{

    struct Point{
        double x;
        double y;

        };

    class Warp {
        public:
        Warp(){}
        virtual ~Warp(){};
        virtual Point operator()(const Point& p) const = 0;

    };

    void ApplyWarp(const Point& cen, const double range, 
        const Warp& w, const ImgProc& in, ImgProc& out);

    class Translate: public Warp{

        public:
            Translate (const Point& t);
            ~Translate(){}

            Point operator()(const Point& P) const;

        private:
            Point trans;
    };
    
 
 
}

