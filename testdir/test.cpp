// CPP program without virtual destructor 
#include<iostream>
 
using namespace std;
 
class base {
public:
  base()     
  { cout<<"Constructing base \n"; }
  void method(){ cout << "base method\n";}
  virtual ~base()
  { cout<<"Destroying base \n"; }     
};
 
class derived: public base {
public:
  derived()     
  { cout<<"Constructing derived \n"; }
  void method(){ cout << "derived method";}
  ~derived()
  { cout<<"Destroying derived \n"; }
};
 
int main(void)
{
  derived *d = new derived();  
  base *b = d;
  delete b;
  return 0;
}
