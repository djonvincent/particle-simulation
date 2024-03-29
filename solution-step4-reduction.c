// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <math.h>
#include <limits>
#include <iomanip>
#include <omp.h>
#include <vector>

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  if (argc == 5) {
    std::string line;
    std::ifstream dataFile (argv[4]);
    if (dataFile.is_open()) {
      getline(dataFile, line);
      dataFile.close();
    } else {
      std::cerr << "Unable to open file " << argv[4] << '\n';
      exit(-2);
    }
    char delim[] = " ";
    char delim2[] = "      ";
    char lineCStr[line.size()+1];
    std::strcpy(lineCStr, line.c_str());
    char * ptr = strtok(lineCStr, delim2);
    std::vector<double*> xVector;
    std::vector<double*> vVector;
    std::vector<double> massVector;
    NumberOfBodies = 0;
    while (ptr != NULL) {
      double* xi = new double[3];
      double* vi = new double[3];
      xi[0] = std::stof(ptr);
      xi[1] = std::stof(strtok(NULL, delim));
      xi[2] = std::stof(strtok(NULL, delim));
      xVector.push_back(xi);

      vi[0] = std::stof(strtok(NULL, delim));
      vi[1] = std::stof(strtok(NULL, delim));
      vi[2] = std::stof(strtok(NULL, delim));
      vVector.push_back(vi);

      massVector.push_back(std::stof(strtok(NULL, delim)));

      if (massVector.back()<=0.0 ) {
        std::cerr << "invalid mass for body " << NumberOfBodies << std::endl;
        exit(-2);
      }
      ptr = strtok(NULL, delim2);
      //printf("%f %f %f %f %f %f %f\n", xi[0], xi[1], xi[2], vi[0], vi[1], vi[2], massVector.back());
      NumberOfBodies ++;
    }
    x    = new double*[NumberOfBodies];
    v    = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];
    for (int i=0; i<NumberOfBodies; i++) {
      x[i] = xVector.at(i);
      v[i] = vVector.at(i);
      mass[i] = massVector.at(i);
    }
  } else {
    NumberOfBodies = (argc-4) / 7;
    x    = new double*[NumberOfBodies];
    v    = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];

    for (int i=0; i<NumberOfBodies; i++) {
      x[i] = new double[3];
      v[i] = new double[3];

      x[i][0] = std::stof(argv[readArgument]); readArgument++;
      x[i][1] = std::stof(argv[readArgument]); readArgument++;
      x[i][2] = std::stof(argv[readArgument]); readArgument++;

      v[i][0] = std::stof(argv[readArgument]); readArgument++;
      v[i][1] = std::stof(argv[readArgument]); readArgument++;
      v[i][2] = std::stof(argv[readArgument]); readArgument++;

      mass[i] = std::stof(argv[readArgument]); readArgument++;

      if (mass[i]<=0.0 ) {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV = 0;
  minDx = std::numeric_limits<double>::max();

  //printf("new timestep\n");
  for (int j=0; j<NumberOfBodies; j++) {
    double f0 = 0;
    double f1 = 0;
    double f2 = 0;
    #pragma omp parallel
    {
      double myf0 = 0;
      double myf1 = 0;
      double myf2 = 0;
      double myMinDx = minDx;
      #pragma omp for
      for (int i=0; i<NumberOfBodies; i++) {
        //printf("j: %d i: %d thread: %d\n", j, i, omp_get_thread_num());
        if (i==j) {
          continue;
        }
        const double distance = sqrt(
          (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
          (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
          (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
        );

        // x,y,z forces acting on particle 0
        double m1m2OverDistanceCubed = mass[i]*mass[j] / (distance * distance * distance);
        myf0 += (x[i][0]-x[j][0]) * m1m2OverDistanceCubed;
        myf1 += (x[i][1]-x[j][1]) * m1m2OverDistanceCubed;
        myf2 += (x[i][2]-x[j][2]) * m1m2OverDistanceCubed;

        myMinDx = std::min( myMinDx,distance );
      }
      #pragma omp critical
      {
        f0 += myf0;
        f1 += myf1;
        f2 += myf2;
        minDx = std::min(minDx, myMinDx);
      }
    }

    v[j][0] += timeStepSize * f0 / mass[j];
    v[j][1] += timeStepSize * f1 / mass[j];
    v[j][2] += timeStepSize * f2 / mass[j];
    maxV = std::max(
      maxV,
      std::sqrt( v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2] )
    );
  }

  for (int i=0; i<NumberOfBodies; i++) {
    x[i][0] += timeStepSize * v[i][0];
    x[i][1] += timeStepSize * v[i][1];
    x[i][2] += timeStepSize * v[i][2];
  }

  // These are three buggy lines of code that we will use in one of the labs
//  x[0][3] = x[0][2] + timeStepSize * v[0][2];
//  x[0][2] = x[0][2] + timeStepSize * v[0][2] / 0.0;
//  x[50000000][1] = x[0][2] + timeStepSize * v[0][2] / 0.0;

  t += timeStepSize;

}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( !((argc-4)%7==0 || argc ==5 )) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
