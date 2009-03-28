#ifndef SHAPELENS_SINGLETON_H
#define SHAPELENS_SINGLETON_H

namespace shapelens {
  /// Meyers Singleton class.
  template <class T>
    class Singleton {
  public:
    /// Get single instance of type T.
    static T& getInstance() {
      static T _instance;
      return _instance;
    }
  private:
    Singleton();          // constructor hidden
    ~Singleton();          // destructor hidden
    Singleton(Singleton const&);    // copy constructor hidden
    Singleton& operator=(Singleton const&);  // assign op hidden
  };

}
#endif
