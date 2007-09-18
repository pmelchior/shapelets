#include <ShapeLens.h>

// convertSIF2PPM:
// create a PPM file from a shapelet model specified by SIF file.

int main(int argc, char *argv[]) {
  if (argc != 6)
    std::cout << "usage: convertSIF2PPM <sif file> <ppm file> <bool: poissonian> <noise mean> <noise sigma>" << std::endl;
  else {
    std::string filename = argv[1];
    std::string ppmname = argv[2];
    int poissonian = atoi(argv[3]);
    data_t noisemean = atof(argv[4]);
    data_t noisesigma = atof(argv[5]);

    // open sif file
    ShapeletObject * sobj =  new ShapeletObject(filename);
    // get grid from sobj
    const Grid& grid = sobj->getGrid();
    // construct shapelet model from coeffs stored in sobj
    NumVector<data_t> data = sobj->getModel();
    // add noise if required
    if (poissonian==1)
      addPoissonianNoise(data,noisemean,noisesigma);
    else
      addGaussianNoise(data,noisemean,noisesigma);

    // write model to ppm file:
    // colormodel in {"RED", "BLUE", "GREEN", "GRAY", "WARM"}
    // data scaling in {"LINEAR", "SQUARE_ROOT", "LOGARITHMIC"} 
    // for the data range between min() and max()
    writePPMImage(ppmname,"SPECTRUM","LINEAR",data.min(),data.max(),grid,data);
    delete sobj;
  }
}
  
							
