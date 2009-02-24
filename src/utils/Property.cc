#include <utils/Property.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

Property::Property() : std::map<std::string, variant_t>() {}

// Helper class for performing write()
class OutFunctor : public boost::static_visitor<void> {
public:
  OutFunctor(std::ostream& o) : out(o) {
  }
  void operator()(int e) const {
    out << "\tI\t" << e << std::endl;
  }
  void operator()(data_t e) const {
    out << "\tD\t" << e << std::endl;
  }
  void operator()(std::string& e) const {
    out << "\tS\t" << e << std::endl;
  }
  void operator()(std::vector<int>& e) const {
    out << "\tVI\t";
    for (std::vector<int>::iterator iter = e.begin(); iter != e.end(); iter++) {
      out << *iter;
      if (iter != --e.end())
	out  << ",";
    }
    out << std::endl;
  }
  void operator()(std::vector<data_t>& e) const {
    out << "\tVD\t";
    for (std::vector<data_t>::iterator iter = e.begin(); iter != e.end(); iter++) {
      out << *iter;
      if (iter != --e.end())
	out  << ",";
    }
    out << std::endl;
  }
private:
  std::ostream & out;
};

void Property::write(std::ostream& o) {
  OutFunctor out(o);
  for (Property::iterator iter = std::map<std::string, variant_t>::begin(); iter != std::map<std::string, variant_t>::end(); iter++) {
    o << iter->first << "\t";
    boost::apply_visitor(out, iter->second);
  }
}

void Property::read(std::istream& in) {
  std::string line, name, type, value;
  typedef boost::tokenizer<boost::char_separator<char> > Tok;
  boost::char_separator<char> fieldsep("\t");
  boost::char_separator<char> vectorsep(",");
  while(getline(in, line)) {
    Tok tok(line, fieldsep);
    int i=0;
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter) {
      switch(i) {
      case 0: name = *tok_iter; break;
      case 1: type = *tok_iter; break;
      case 2: value = *tok_iter; break; 
      }
      i++;
    } 
    if (i==3) { // name, type, value filled
      if (type[0] == 'V') { // vector type
	Tok vtok(value, vectorsep);
	if (type[1] == 'I') {
	  std::vector<int> vi;
	  for(Tok::iterator tok_iter = vtok.begin(); tok_iter != vtok.end(); ++tok_iter) 
	    vi.push_back(boost::lexical_cast<int>(*tok_iter));
	  std::map<std::string, variant_t>::insert(std::make_pair(name,vi));
	} else if (type[1] == 'D') {
	  std::vector<data_t> vi;
	  for(Tok::iterator tok_iter = vtok.begin(); tok_iter != vtok.end(); ++tok_iter) 
	    vi.push_back(boost::lexical_cast<data_t>(*tok_iter));
	  std::map<std::string, variant_t>::insert(std::make_pair(name,vi));
	}
      } else {
	if (type[0] == 'S')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,value));
	else if (type[0] == 'I')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,boost::lexical_cast<int>(value)));
	else if (type[0] == 'D')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,boost::lexical_cast<data_t>(value)));
      }
    }
  }
}
