

#include <iomanip>
#include <complex>
#include <cmath>
#include "warp.h"
using namespace std;

namespace img {
    class ColorLUT{

        public:

            ColorLUT(double gamma = 1.0);

            ~ColorLUT(){};

            void operator()(const double& value, vector<float>& C) const;
        private:
            double gamma;
            vector<float>  black;
            vector<vector<float>> bands;

    };

    void ApplyFractalWarp(const Point& center, const double range, const Warp& w, const ColorLUT& lut, ImgProc& out);

 
    class JuliaSet : public Warp
    {
    public:
        JuliaSet(const Point &PC, const int nb, const int cycle);
        ~JuliaSet(){};

        Point operator()(const Point &P) const;

    private:
        Point c;
        int nb_interations;
        int cycles;
    };

    void do_cylinder_IFS(ImgProc &out, vector<Point> &p_vec, ColorLUT &lut);
    // namespace img
    void do_linear_IFS(ColorLUT &lut, ImgProc &out, vector<Point> &p_vec);
    void linear_IFS(int ic, const Point &p, vector<Point> &vec_p, ColorLUT &lut, ImgProc &out);
    void power_IFS(int ic, const Point &p, vector<Point> &vec_p, ColorLUT& lut, ImgProc &out, double f_theta);
    void do_power_IFS(ImgProc &out, vector<Point> &p_vec, ColorLUT& lut, double fractal_theta);
    void cylinder_IFS(const Point &p, vector<Point> &p_vec, ColorLUT &lut, ImgProc &out);
}
