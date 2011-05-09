#include <shapelens/ShapeLens.h>
#include <shapelens/shapelets/ShapeletObject.h>
#include <tclap/CmdLine.h>

using namespace shapelens;
using namespace std;

int main(int argc, char *argv[]) {
  TCLAP::CmdLine cmd("Convolve ShapeletObject with anoher ShapeletObject", ' ', "0.3");
  TCLAP::SwitchArg deconv("d","deconvolve","Deconvolve from kernel", cmd, false);
  TCLAP::ValueArg<string> input("i","input","input SIF file",true,"","string", cmd);
  TCLAP::ValueArg<string> kernel ("k","kernel","kernel SIF file",true,"","string",cmd);
  TCLAP::ValueArg<string> output ("o","output","output SIF file",true,"","string",cmd);
  TCLAP::ValueArg<unsigned int> nmax_f ("n","nmax_f","n_max of the deconvolved object",false,0,"unsigned int",cmd);
  cmd.parse(argc,argv);
  
  // open kernel file and normalize
  ShapeletObject sk(kernel.getValue());
  sk.brighten(1./sk.integrate());

  ShapeletObject si(input.getValue());

  if (deconv.isSet()) {
    if (nmax_f.isSet()) {
      ImageTransformation trafo;
      CoefficientVector<data_t> coeffs = si.getCoeffs();
      data_t beta_f = sqrt(si.getBeta()*si.getBeta() - sk.getBeta()*sk.getBeta());
      NumMatrix<data_t> P = trafo.getConvolutionMatrix(sk.getCoeffs(),nmax_f.getValue(), sk.getNMax(),si.getNMax(),beta_f,sk.getBeta(),si.getBeta());
      coeffs = (P.transpose()*P).invert()*(P.transpose()*coeffs);
      si.setCoeffs(coeffs);
      si.setBeta(beta_f);
    } 
    else
      si.deconvolve(sk.getCoeffs(),sk.getBeta());
  }
  else
    si.convolve(sk.getCoeffs(),sk.getBeta());

  si.save(output.getValue());
}
