#include <frame/Profile.h>
#include <gsl/gsl_math.h>

Profile::Profile() {
  pixels = 0;
}

//Profile::Profile(Point2D& center) {
//}
  
Profile::Profile(const Point2D& instart, const Point2D& instop) {
  start = instart;
  stop = instop;
  center = Point2D((stop(0)-start(0))/2 + start(0)-1,(stop(1)-start(1))/2 + start(1)-1); 
  pixels = GSL_MAX_INT((int)(floor(stop(0)) - floor(start(0))), (int) (floor(stop(1)) - floor(start(1))));
  distance_value.resize(pixels,2);
}

data_t Profile::getDistance(unsigned int i) {
  return distance_value(i,0);
}

data_t Profile::getValue(unsigned int i) {
  return distance_value(i,1);
}

Point2D& Profile::getStart() {
  return start;
}

Point2D& Profile::getEnd() {
  return stop;
}

Point2D& Profile::getCenter() {
  return center;
}

int Profile::size() {
  return pixels;
}

void Profile::calculate(NumVector<data_t>& data, int axsize) {
  int diffx = (int) (floor(stop(0)) - floor(start(0)));
  int diffy = (int) (floor(stop(1)) - floor(start(1)));
  if (diffx == diffy) {
    for (int i = 0; i < pixels; i++) {
      int x = (int)floor(start(0)) + i;
      int y = (int)floor(start(1)) + i;
      ///std::cout << x  <<  " " << y << std::endl;
      unsigned int index = y*axsize + x;
      data_t distance = sqrt((x - floor(center(0)))*(x - floor(center(0))) +(y - floor(center(1)))*(y - floor(center(1))));
      if (x < floor(center(0)) || y < floor(center(1))) distance *= -1;
      distance_value(i,0) = distance;
      distance_value(i,1) = data(index);
    }
  }
  if (diffx > diffy) {
    for (int i = 0; i < pixels; i++) {
      int x = (int)floor(start(0)) + i;
      int y = (int)floor(start(1)) + i*diffy/diffx;
      unsigned int index = y*axsize + x;
      data_t distance = sqrt((x - floor(center(0)))*(x - floor(center(0))) +(y - floor(center(1)))*(y - floor(center(1))));
      if (x < floor(center(0)) || y < floor(center(1))) distance *= -1;
      distance_value(i,0) = distance;
      distance_value(i,1) = data(index);
    }
  }
  if (diffy > diffx) {
    for (int i = 0; i < pixels; i++) {
      int x = (int)floor(start(0)) + i*diffx/diffy;
      int y = (int)floor(start(1)) + i;
      unsigned int index = y*axsize + x;
      data_t distance = sqrt((x - floor(center(0)))*(x - floor(center(0))) +(y - floor(center(1)))*(y - floor(center(1))));
      if (x < floor(center(0)) || y < floor(center(1))) distance *= -1;
      distance_value(i,0) = distance;
      distance_value(i,1) = data(index);
    }
  }
}
