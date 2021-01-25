//------------------------------------------------
//
//  img_paint
//
//
//-------------------------------------------------




#include <cmath>
#include <omp.h>
#include "imgproc.h"
#include "CmdLineFind.h"
#include "Color_LUT.h"

#include <vector>
#include <fstream>
#include <iomanip>

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include <iostream>
#include <stack>


using namespace std;
using namespace img;

ImgProc image;
vector<ImgProc> undo_stack;
vector<unsigned char> key_stack;


//double fractal_theta = 0.0;
int julia_iterations = 10;
int julia_cycles = 10;
double fractal_theta = 0.6;
vector<float> min_c;
vector<float> max_c;
vector<float> mean_c;
vector<float> standard_dev_c;
vector<vector<float>> histg1;

vector<vector<float>> CDF1;
vector<vector<int>> histg2;
vector<vector<float>> his_pix;

vector<vector<float>> CDF2;
vector<vector<float>> PDF1;
vector<vector<float>> PDF2;

ofstream data;

//float


void setNbCores( int nb )
{
   omp_set_num_threads( nb );
}

void cbMotion( int x, int y )
{
}

void cbMouse( int button, int state, int x, int y )
{
}

void cbDisplay( void )
{
   glClear(GL_COLOR_BUFFER_BIT );
   glDrawPixels( image.nx(), image.ny(), GL_RGB, GL_FLOAT, image.raw() );
   glutSwapBuffers();
}

void cbIdle()
{
   glutPostRedisplay();
}


void deFractal(const Warp& s) {
  ColorLUT lut;

  Point center;
  center.x = 0.3001*0.8;
  center.y = 0.1505*0.8;
  double range = 2;
  img::ApplyFractalWarp(center, range, s, lut, image);
  glutPostRedisplay();

}



void cbOnKeyboard( unsigned char key, int x, int y )
{


  int julia_iterations = 100;
  int julia_cycles = 30;
  //double fractal_theta = 0.0;
  switch (key)
  {

  case 'i': {

    cout << "linear IFS\n";
    ColorLUT lut;
    vector<Point> p_vec;

    img::do_linear_IFS(lut, image, p_vec);

    }
      break;

    case 'c' : 
    image.compliment();
    cout << "Compliment\n";

    break;
  case 'V':
    //brit *=1.1;
    image.brightness(1.05);
    cout << "brightness increase\n";
    break;
  case 'v':
    //  brit *= 0.90;
    image.brightness(0.95);
    cout << "brightness decrease\n";
    break;
  case 'B':

    image.bias(0.05);
    cout << "bias increase\n";
    break;
  case 'b':

    image.bias(-0.05);
    cout << "bias decrease\n";
    break;

  case 'f':
    image.flap();
    cout << "flap\n";
    break;

  case 'G':
    // gm*=1.1;
    image.gamma(1.8);

    cout << "increase gamma\n";
    break;
  case 'g':

    image.gamma(0.95);
    cout << "decrease gamma\n";
    break;
  case 'q':
    image.quantize();
    cout << "quantize\n";
    break;
  case 'w':
    image.grayscale();
    cout << "grayscale\n";
    break;
  case 'C':
    image.rms();
    cout << "rms contrast\n";
    break;
  case 'r':
    image.org_data();
    cout << "original data\n"
         << endl;
    break;
  case 'o':
    image.write();
    cout << "write file\n";
    break;
  case 'd':
    data.open("analysis.txt");
    //image.gamma(1.8);

    data << "\n\n		PDF\n\n     Red    	   green     	 Blue         after his_eq    		Red          Green    	Blue\n";
    image.histogram(histg1, CDF1, PDF1, his_pix);
    his_pix.clear();
    histg1.clear();
    image.histogram(histg1, CDF2, PDF2, his_pix);
    //original RBG
    for (long i = 0; i < 255; i++)
    {

      for (long c = 0; c < 3; c++)
      {

        data << "   " << setw(9) << PDF1[c][i] << "    ";
      }

      data << "     ||     ";
      for (long c = 0; c < 3; c++)
      {

        data << "      " << setw(9) << PDF2[c][i];
      }
      data << endl;
    }
    cout << endl;

    data << endl;

    CDF1.clear();

    PDF1.clear();
    CDF2.clear();

    PDF2.clear();
    his_pix.clear();
    histg1.clear();

    data << "\n\n\n-------------------------------------------with 1.8 gamma --------------------------\n\n\n"
         << endl;
    data << "			PDF\n     Red    	  green  	 Blue     	     after his_eq		Red       	Green   	Blue\n";

    image.org_data();
    image.flap();

    his_pix.clear();
    //CDF2.clear();
    //PDF2.clear();
    image.gamma(1.8);
    image.histogram(histg1, CDF1, PDF1, his_pix);
    his_pix.clear();
    image.histogram(histg1, CDF2, PDF2, his_pix);

    //original RBG
    for (long i = 0; i < 255; i++)
    {

      for (long c = 0; c < 3; c++)
      {

        data << "   " << setw(9) << PDF1[c][i] << "    ";
      }

      data << "       ||    ";
      for (long c = 0; c < 3; c++)
      {

        data << "      " << setw(9) << PDF2[c][i];
      }
      data << endl;
    }

    data.close();
    his_pix.clear();
    histg1.clear();
    //histg2.clear();
    CDF1.clear();
    //CDF2.clear();
    PDF1.clear();
    break;

  case 'h':
    image.histogram(histg1, CDF1, PDF1, his_pix);
    histg1.clear();
    PDF1.clear();

    CDF1.clear();
    cout << "histogram\n";

    break;
  case 's':

    cout << "stat \n";

    image.stat();
    break;

  case 'J' :{
      cout << "julia\n";

      Point JuliaP;
      //----------- V13(x, y) = √r · (cos(θ/2 + Ω),sin(θ/2 + Ω)
      undo_stack.push_back(image);
      key_stack.push_back(key);
      fractal_theta += 1.0;
      JuliaP.x = 0.8 * cos(fractal_theta + 3.14159 / 100.0);
      JuliaP.y = 0.1 * sin(fractal_theta + 3.14159 / 100.0);
      deFractal(JuliaSet(JuliaP, julia_iterations, julia_cycles));
  }
    break;


  case 'p': 
    {
      cout << "power IFS\n";
      ColorLUT lut;
      vector<Point> p_vec;

      //img::o_power_IFS(lut, image, p_vec);
      fractal_theta +=0.5;
      do_power_IFS(image, p_vec, lut, fractal_theta);
    }
    break;
  case 'z' : {

    cout << "cylinder IFS\n";
    //V29(x, y) = (sin x, y)
    ColorLUT lut;
    vector<Point> p_vec;
    do_cylinder_IFS(image, p_vec, lut);
   
  }
  break;
 }
}

void PrintUsage()
{
   cout << "img_paint keyboard choices\n";
   cout << "c         compliment\n";
   cout << "V         brightness increase\n";
   cout << "v         brightness decrease\n";
   cout << "B         bias increase\n";
   cout << "b         bias decrease\n";
   cout << "f         flap\n";
   cout << "G         increase gamma\n";
   cout << "g         decrease gamma\n";
   cout << "q         quantize\n";
   cout << "w         grayscale\n";
   cout << "R         original data\n";
   cout << "o         write out_file \n" << endl;
   cout << "h         histogram\n";
   cout << "s	        statistics\n";
   cout << "d	        analysis.txt\n";
   cout << "i         linear IFS\n";
   cout << "p         power  IFS\n";
   cout << "J         JuliaSet\n";
   cout << "z         cylinder\n";

}


int main(int argc, char** argv)
{
   lux::CmdLineFind clf( argc, argv );

   setNbCores(8);

   string imagename = clf.find("-image", "", "Image to drive color");

   clf.usage("-h");
   clf.printFinds();
   PrintUsage();

   image.load(imagename);


   // GLUT routines
   glutInit(&argc, argv);

   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   glutInitWindowSize( image.nx(), image.ny() );

   // Open a window
   char title[] = "img_paint";
   glutCreateWindow( title );

   glClearColor( 1,1,1,1 );

   glutDisplayFunc(&cbDisplay);
   glutIdleFunc(&cbIdle);
   glutKeyboardFunc(&cbOnKeyboard);
   glutMouseFunc( &cbMouse );
   glutMotionFunc( &cbMotion );

   glutMainLoop();
   return 1;
};
