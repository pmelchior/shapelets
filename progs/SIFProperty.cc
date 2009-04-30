#include <ShapeLens.h>
#include <tclap/CmdLine.h>
#include <fstream>

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Read/Modify ShapeletObject's property struct", ' ', "0.1");
  TCLAP::ValueArg<std::string> input("i","input","File with Property data",false,"","string",cmd);
  TCLAP::ValueArg<std::string> output("o","output","File to write Property data to",false,"","string",cmd);
  TCLAP::SwitchArg clear("c","clear","Clear existing property from SIF file", cmd, false);
  TCLAP::UnlabeledValueArg<std::string> siffile("file","SIF file to analyze",true,"","string", cmd);
  cmd.parse(argc,argv);

  shapelens::ShapeletObject sobj(siffile.getValue());
  if (clear.isSet()) {
    sobj.prop.clear();
    if (!input.isSet())
      sobj.save(siffile.getValue());
  }
  if (input.isSet()) {
    std::ifstream ifs(input.getValue().c_str());
    sobj.prop.read(ifs);
    ifs.close();
    sobj.save(siffile.getValue());
  }
  if (output.isSet()) {
    std::ofstream ofs(output.getValue().c_str());
    sobj.prop.write(ofs);
    ofs.close();
  } else
    sobj.prop.write(std::cout);
}
