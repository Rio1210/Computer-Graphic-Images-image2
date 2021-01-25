#include "Color_LUT.h"
#include <iostream>

using namespace img;
using namespace std;

ColorLUT::ColorLUT(double gam) : gamma(gam)
{

    vector<float> C;
    C.push_back(0.0);
    C.push_back(0.0);
    C.push_back(0.0);
    black = C;
    //yuyan yingyuanse
    C[0] = 150.0 / 255.0;
    C[1] = 230.0 / 255.0;
    C[2] = 161.0 / 255.0;
    bands.push_back(C);
    //show yingyuanse
    C[0] = 154.0 / 255.0;
    C[1] = 219.0 / 255.0;
    C[2] = 232.0 / 255.0;
    bands.push_back(C);
    //blue
    C[0] = 58.0 / 255.0;
    C[1] = 73.0 / 255.0;
    C[2] = 88.0 / 255.0;
    bands.push_back(C);
    //tillman brick
    C[0] = 86.0 / 255.0;
    C[1] = 46.0 / 255.0;
    C[2] = 25.0 / 255.0;
    bands.push_back(C);

    //game day sky
    C[0] = 156.0 / 255.0;
    C[1] = 157.0 / 255.0;
    C[2] = 192.0 / 255.0;
    bands.push_back(C);

    //pink
    C[0] = 250.0 / 255.0;
    C[1] = 182.0 / 255.0;
    C[2] = 193.0 / 255.0;
    bands.push_back(C);

    //gray
    C[0] = 219.0 / 255.0;
    C[1] = 219.0 / 255.0;
    C[2] = 219.0 / 255.0;
    bands.push_back(C);

    //arbi
    C[0] = 100.0 / 255.0;
    C[1] = 100.0 / 255.0;
    C[2] = 100.0 / 255.0;
    bands.push_back(C);
}

void ColorLUT::operator()(const double &value, vector<float> &C) const
{

    C = black;

    if (value > 1.0 || value < 0.0)
    {
        return;
    }

    double x = pow(value, gamma) * (bands.size() - 1);

    size_t low_index = (size_t)x;
    size_t high_index = low_index + 1;
    double weight = x - (double)low_index;

    if (high_index >= bands.size())
    {
        high_index = bands.size() - 1;
    }

    for (size_t c = 0; c < C.size(); c++)
    {
        C[c] = bands[low_index][c] * (1.0 - weight) 
            + bands[high_index][c] * weight;
    }
}

void img::do_power_IFS(ImgProc &out, vector<Point> &p_vec, ColorLUT &lut, double fractal_theta){

    Point P;

    P.x = 2.0 * drand48() - 1; //(-1 ~ 1)
    P.y = 2.0 * drand48() - 1;
    float w = drand48();
    // cout << "check\n";
    for (size_t iter = 0; iter < 50; iter++)
    {
        size_t ic = (size_t)(drand48() * 100);

        power_IFS(ic, P, p_vec, lut,  out, fractal_theta);
        P = p_vec[ic];
        w = (w * out.img_data[out.index(P.x, P.y, 0)]) / 2; //p_vec[ic])
      
        if (iter > 20)
        {
 
            if (P.x >= -1.0 && P.x <= 1.0 && P.y >= -1.0 && P.y <= 1.0)
            {
                int ii = (int)P.x * out.nx();

                int jj = (int)P.y * out.ny();

                if (ii < out.nx())
                {

                    if (jj < out.ny())
                    {
                        vector<float> color;
                        lut(w, color);
                        vector<float> cc(3);
                        out.value(jj, ii, cc);
                        #pragma omp parallel for
                        for (size_t iic = 0; iic < cc.size() - 1; iic++)
                        {
                            //cout << "in loop loop loop\n";                      
                            cc[iic] *= cc[cc.size() - 1];
                            cc[iic] += color[iic];
                            cc[iic] /= cc[cc.size() - 1];
                            //cc[iic] = (cc[iic] + color[iic]) / (cc[cc.size() - 1] + 1);
                        }
                        cc[cc.size() - 1] += 1;
                        out.set_value(ii, jj, cc);
                    }
                }
            }
        }
    }
}


void img::power_IFS(int ic, const Point &p, vector<Point> &vec_p, ColorLUT& lut,  ImgProc &out, double fractal_theta)
{

    Point pp;
    //V19(x, y) = pow(r,sin θ) * (cos θ, sin θ)

    for (int j = 0; j < out.ny(); j++)
    {
        //#pragma omp parallel for
        for (int i = 0; i < out.nx(); i++)
        {

            double rate = sqrt((pp.x * pp.x) + (pp.y * pp.y)) / 2;

            pp.x += pow(rate, sin(fractal_theta))*(cos(fractal_theta));
            pp.y += pow(rate, sin(fractal_theta)) * (sin(fractal_theta));


            vector<float> v(3, 0.0);

            lut(rate, v);
            out.set_value(i, j, v);

            vec_p.push_back(pp);
            // cout << "check 5 pushback loop inlear\n";
        }
    }
}

//----------- V13(x, y) = √r · (cos(θ/2 + Ω),sin(θ/2 + Ω)
void img::ApplyFractalWarp(const Point &center, const double range,
                           const Warp &w, const ColorLUT &lut, ImgProc &out)
{

    float R = 2.0;

    for (int j = 0; j < out.ny(); j++)
    {
#pragma omp parallel for
        for (int i = 0; i < out.nx(); i++)
        {
            Point P;
            P.x = 2.0 * (double)i / (double)out.nx() - 1.0;
            P.y = 2.0 * (double)j / (double)out.ny() - 1.0;
            P.x *= range;
            P.y *= range;
            P.x += center.x;
            P.y += center.y;
            Point PP = w(P);

            double rate = sqrt((PP.x * PP.x) + (PP.y * PP.y)) / R;
            vector<float> v(3, 0.0);
            lut(rate, v);
            out.set_value(i, j, v);
            //cout << "in apply P.x:  " << P.x << endl;
        }
    }
}

//----------

JuliaSet::JuliaSet(const Point& P0, const int nb, const int cycle) 
    : c(P0), nb_interations(nb), cycles(cycle){}

Point JuliaSet::operator() ( const Point &P ) const {
    std::complex<double> Pc(P.x, P.y);
    std::complex<double> CC(c.x, c.y);

    for(int i = 0; i < nb_interations; i++) {
        std::complex<double> temp = Pc;
        for(int c = 1; c < cycles; c++) {
            temp  = temp * Pc;

        }
        Pc = temp + CC;
    }

    Point pout;
    pout.x = Pc.real();
    pout.y = Pc.imag();
    return pout;
}

void img::linear_IFS(int ic, const Point &p, vector<Point> &vec_p, ColorLUT &lut, ImgProc &out)
{

    Point pp;
   
    for (int j = 0; j < out.ny(); j++)
    {
//#pragma omp parallel for
        for (int i = 0; i < out.nx(); i++)
        {
            pp.x += p.x;// / (double)out.nx();
            pp.y += p.y;// / (double) out.ny();
     
            vector<float> v(3, 0.0);
         //   double rate = sqrt((pp.x * pp.x) + (pp.y * pp.y)) / 2;
           // lut(rate, v);
            out.set_value(i, j, v);
        //
            vec_p.push_back(pp);
           // cout << "check 5 pushback loop inlear\n";
        }
    }
}

void img::do_linear_IFS(ColorLUT &lut, ImgProc &out, vector<Point> &p_vec)
{
    Point P;
  //  Warp warp;
    //ColorLUT lut;
   // ImgProc &out 
    P.x = 2.0 * drand48() - 1; //(-1 ~ 1)
    P.y = 2.0 * drand48() - 1;
    // float w = drand48()*0.5;//+10.0)/2;
    float w = drand48();
   // cout << "check\n";
    for (size_t iter = 0; iter < 50; iter++)
    {
        size_t ic = (size_t)(drand48() * 100);
       

        linear_IFS(ic, P, p_vec, lut, out);
        P = p_vec[ic];

        w = (w * out.img_data[out.index(P.x, P.y, 0)]) / 2; //p_vec[ic])

        if (iter > 20)
        {

      
            if (P.x >= -1.0 && P.x <= 1.0 && P.y >= -1.0 && P.y <= 1.0)
            {
                int ii = (int)P.x * out.nx();

                int jj = (int)P.y * out.ny();

                if (ii < out.nx())
                {

                    if (jj < out.ny())
                    {
                        vector<float> color;
                        lut(w, color);
                        vector<float> cc(3);
                        out.value(jj, ii, cc);
                        for (size_t iic = 0; iic < cc.size() - 1; iic++)
                        {
                            ;
                            cc[iic] *= cc[cc.size() - 1];
                            cc[iic] += color[iic];
                            cc[iic] /= cc[cc.size() - 1];
                            //cc[iic] = (cc[iic] + color[iic]) / (cc[cc.size() - 1] + 1);
                        }
                        cc[cc.size() - 1] += 1;
                        out.set_value(ii, jj, cc);
                    }
                }
            }
            
        }
    }
}

void img::do_cylinder_IFS(ImgProc &out, vector<Point> &p_vec, ColorLUT &lut){

    Point P;
    //  Warp warp;
    //ColorLUT lut;
    // ImgProc &out
    P.x = 2.0 * drand48() - 1; //(-1 ~ 1)
    P.y = 2.0 * drand48() - 1;
    // float w = drand48()*0.5;//+10.0)/2;
    float w = drand48();
    // cout << "check\n";
    for (size_t iter = 0; iter < 100; iter++)
    {
        size_t ic = (size_t)(drand48() * 100);

        cylinder_IFS(P, p_vec, lut,out);
        P = p_vec[ic];

        w = (w * out.img_data[out.index(P.x, P.y, 0)]) / 2; //p_vec[ic])

        if (iter > 20)
        {

            if (P.x >= -1.0 && P.x <= 1.0 && P.y >= -1.0 && P.y <= 1.0)
            {
                int ii = (int)P.x * out.nx();

                int jj = (int)P.y * out.ny();

                if (ii < out.nx())
                {

                    if (jj < out.ny())
                    {
                        vector<float> color;
                        lut(w, color);
                        vector<float> cc(3);
                        out.value(jj, ii, cc);
                        for (size_t iic = 0; iic < cc.size() - 1; iic++)
                        {
                            
                            cc[iic] *= cc[cc.size() - 1];
                            cc[iic] += color[iic];
                            cc[iic] /= cc[cc.size() - 1];
                            //cc[iic] = (cc[iic] + color[iic]) / (cc[cc.size() - 1] + 1);
                        }
                        cc[cc.size() - 1] += 1;
                        out.set_value(ii, jj, cc);
                    }
                }
            }
        }
    }
}

void img::cylinder_IFS(const Point &p, vector<Point> &p_vec, ColorLUT &lut, ImgProc &out)
{

    Point pp;

    for (int j = 0; j < out.ny(); j++)
    {
        //#pragma omp parallel for
        for (int i = 0; i < out.nx(); i++)
        {
            pp.x += sin(p.x); // / (double)out.nx();
            pp.y += p.y; // / (double) out.ny();

            vector<float> v(3, 0.0);
            //   double rate = sqrt((pp.x * pp.x) + (pp.y * pp.y)) / 2;
            //  lut(rate, v);
            out.set_value(i, j, v);

            p_vec.push_back(pp);
            // cout << "check 5 pushback loop inlear\n";
        }
    }
}