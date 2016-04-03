#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
using namespace std;

class Myclass
{
public:
	double r;
	Myclass()
	{r=0; s=0;} 
	Myclass(double mr, double ms);
	~Myclass() {}
	double add();
	double minus();
private:
	double s;
};

Myclass::Myclass(double mr, double ms)
{
	r=mr;
	s=ms;
}

double Myclass::add()
{	
	double result;
	result=r+s;
	cout<<result<<endl;
}

double Myclass::minus()
{
	double result;
	result=r-s;
 	cout<<result<<endl;
}

int main(void)
{
	Myclass x(10.0, 20.0);
	x.add();	
	x.minus();
}
