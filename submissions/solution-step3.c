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
#include <math.h>
#include <limits>
#include <iomanip>


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
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

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
  // Exit early and output position of last body if 1 body is left
  if (NumberOfBodies == 1) {
    std::cout << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
    closeParaviewVideoFile();
    exit(0);
  }
  const int numBuckets = 10;
  const double vBucket = maxV/numBuckets;
  int buckets[NumberOfBodies];

  // Sort the objects into 10 buckets based on current velocity such that the
  // fastest object is in the 10th bucket
  for (int i=0; i<NumberOfBodies; i++) {
    if (vBucket > 0) {
      const double vSize = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      buckets[i] = (int)(vSize/vBucket);
      if (vSize >= maxV) {
        buckets[i] = numBuckets-1;
      }
    } else {
      buckets[i] = 0;
    }
  }

  maxV   = 0.0;
  double maxVSquared = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // New positions of bodies
  double newx0[NumberOfBodies];
  double newx1[NumberOfBodies];
  double newx2[NumberOfBodies];

  for (int bucket=0; bucket<numBuckets; bucket++) {
    // Run 2^(bucket no.) timesteps on each bucket
    const int timeSteps = pow(2, bucket);
    const double deltaT = timeStepSize/timeSteps;
    for (int ts=0; ts<timeSteps; ts++) {
      for (int j=0; j<NumberOfBodies; j++) {
        if (buckets[j] != bucket) {
          continue;
        }
        double f0 = 0;
        double f1 = 0;
        double f2 = 0;
        for (int i=0; i<NumberOfBodies; i++) {
          if (i == j) {
            continue;
          }
          const double distance = sqrt(
            (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
            (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
            (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
          );

          // Calculate the forces between object i and j
          // Add forces to object j
          double m1m2OverDistanceCubed = mass[i]*mass[j] / (distance * distance * distance);
          f0 += (x[i][0]-x[j][0]) * m1m2OverDistanceCubed;
          f1 += (x[i][1]-x[j][1]) * m1m2OverDistanceCubed;
          f2 += (x[i][2]-x[j][2]) * m1m2OverDistanceCubed;

          // Only update the minDx for the last timestep
          if (ts == timeSteps-1) {
            minDx = std::min( minDx,distance );
          }
        }

        // Update velocity of object j
        v[j][0] += deltaT * f0 / mass[j];
        v[j][1] += deltaT * f1 / mass[j];
        v[j][2] += deltaT * f2 / mass[j];

        // Update new position of object j
        newx0[j] = x[j][0] + deltaT * v[j][0];
        newx1[j] = x[j][1] + deltaT * v[j][1];
        newx2[j] = x[j][2] + deltaT * v[j][2];
      }

      // Perform continuous collision detection on bodies in the current bucket
      // against all other bodies 
      int collisions[NumberOfBodies]; // Stores bodies that collide
      for (int i=0; i<NumberOfBodies; i++) {
        collisions[i] = -1; // -1 indicates no collision
      }
      double tCollides[NumberOfBodies]; // Stores the times of collisions
      for (int i=0; i<NumberOfBodies; i++) {
        for (int j=i+1; j<NumberOfBodies; j++) {
          // Skip if either object is not in the current bucket
          if (buckets[i] != bucket && buckets[j] != bucket) {
            continue;
          }
          // Solve for the time at which the distance between object i and j is
          // twice the radius (1e-2)
          const double a = (v[i][0]-v[j][0]) * (v[i][0]-v[j][0])  + 
            (v[i][1]-v[j][1]) * (v[i][1]-v[j][1]) +
            (v[i][2]-v[j][2]) * (v[i][2]-v[j][2]);
          const double b = 2*(
            (x[i][0]-x[j][0]) * (v[i][0]-v[j][0]) +
            (x[i][1]-x[j][1]) * (v[i][1]-v[j][1]) +
            (x[i][2]-x[j][2]) * (v[i][2]-v[j][2])
          );
          const double c = (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
            (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
            (x[i][2]-x[j][2]) * (x[i][2]-x[j][2]) -
            (1e-2)*(1e-2);
          const double det = b*b - 4*a*c;
          if (det < 0) { // No solution
            continue;
          }
          double sqrtDet = sqrt(det);
          double tCollide = (-b-sqrtDet)/(2*a);
          if (tCollide < 0) { // If time < 0, use the other solution
            tCollide = (-b+sqrtDet)/(2*a);
          }
          // If 0 <= t <= time step size, then collision happens
          if (tCollide >= 0 && tCollide <= deltaT) {
            collisions[i] = j;
            tCollides[i] = tCollide;
            break;
          }
        }
      }

      // Caculate positions and velocities of fused particles
      for (int i=0; i<NumberOfBodies; i++) {
        int j = collisions[i];
        if (j == -1) {
          continue;
        }
        double tCollide = tCollides[i];
        const double frac = mass[j] / (mass[i]+mass[j]);
        v[i][0] = frac * v[j][0] + (1-frac) * v[i][0];
        v[i][1] = frac * v[j][1] + (1-frac) * v[i][1];
        v[i][2] = frac * v[j][2] + (1-frac) * v[i][2];
        mass[i] = mass[i] + mass[j];
        // Use the contact point of the collision as the new position
        newx0[i] = (x[j][0] + x[i][0] + (v[j][0] + v[i][0])*tCollide) / 2;
        newx1[i] = (x[j][1] + x[i][1] + (v[j][1] + v[i][1])*tCollide) / 2;
        newx2[i] = (x[j][2] + x[i][2] + (v[j][2] + v[i][2])*tCollide) / 2;
        // Move object into the current bucket so that it receives the correct
        // number of updates
        buckets[i] = bucket;

        // Remove object j
        NumberOfBodies--;
        for (int k=j; k<NumberOfBodies; k++) {
          x[k][0] = x[k+1][0];
          x[k][1] = x[k+1][1];
          x[k][2] = x[k+1][2];
          newx0[k] = newx0[k+1];
          newx1[k] = newx1[k+1];
          newx2[k] = newx2[k+1];
          v[k][0] = v[k+1][0];
          v[k][1] = v[k+1][1];
          v[k][2] = v[k+1][2];
          mass[k] = mass[k+1];
          buckets[k] = buckets[k+1];
          collisions[k] = collisions[k+1];
          tCollides[k] = tCollides[k+1];
        }
      }

      // Update the positions and maxV of all particles in the current bucket
      for (int i=0; i<NumberOfBodies; i++) {
        if (buckets[i] == bucket) {
          x[i][0] = newx0[i];
          x[i][1] = newx1[i];
          x[i][2] = newx2[i];
          // Only update the maxV for the last timestep
          if (ts == timeSteps-1) {
            maxVSquared = std::max(
              maxVSquared,
              v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]
            );
          }
        }
      }
    }
  }

  maxV = std::sqrt(maxVSquared);
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
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
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
