#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

class atom 
{
public:
	int atomIndex;
	string atomName;
	string atomType;
	double atomMass;
	double atomCharge;
	double atomUnset;
	atom()
	{
	atomIndex=0;
	atomName="";
	atomType="";
	atomMass=0.0;
	atomCharge=0.0;
	atomUnset=0.0;
	}
	void readAtomLine(const string &text);
	int atomOutputIndex(int shift);
	const string atomOutputName();
	const string atomOutputType();
	double atomOutputMass();
	double atomOutputCharge();
	double atomOutputUnset();
	void atomOutputAtom(int shift);
	~atom(){};
private:
};

void atom::readAtomLine(const string &text)
{
	string atomTitle;
	//sscanf(text, "%s%d%s%s%lf%lf%lf", &atomTitle, &atomIndex, &atomName, &atomType, &atomMass, &atomCharge, &atomUnset);
	istringstream str(text);
	str >> atomTitle >> atomIndex >> atomName >> atomType >> atomMass >> atomCharge >> atomUnset;

}

int atom::atomOutputIndex(int shift)
{
	//printf("%10d\n", atomIndex+shift);
	return(atomIndex+shift);
	//cout<<atomIndex+shift<<endl;
}

const string atom::atomOutputName()
{
	//printf("%10s\n", atomName);
	return(atomName);
}

const string atom::atomOutputType()
{
	//printf("%10s\n", atomType);
	return(atomType);
}

double atom::atomOutputMass()
{
	//printf("%10.4f\n", atomMass);
	return(atomMass);
}

double atom::atomOutputCharge()
{
	//printf("%10.4f\n", atomCharge);
	return(atomCharge);
}

double atom::atomOutputUnset()
{
	//printf("%10.4f\n", atomUnset);
	return(atomUnset);
}

void atom::atomOutputAtom(int shift)
{
	/*printf("%10d\t%10s\t%10s\t%10.4f\t%10.4f\t%10.4f\n",
		atomOutputIndex(shift), atomOutputName(), atomOutputType(),
		atomOutputMass(), atomOutputCharge(), atomOutputUnset());	*/
	cout<<setw(10)<<setprecision(10)<<atomOutputIndex(shift);
	cout<<setw(10)<<atomOutputName();
	cout<<setw(10)<<atomOutputType();
	cout<<setw(10)<<setprecision(4)<<setiosflags(ios::fixed)<<atomOutputMass();
	cout<<setw(10)<<setprecision(2)<<setiosflags(ios::fixed)<<atomOutputCharge();
	cout<<setw(10)<<setprecision(1)<<setiosflags(ios::fixed)<<atomOutputUnset()<<endl;
}

/*int main(void)
{
	vector<string> lines;
	lines.push_back(string("ATOM       1      OA          OA        31.0337    0.00   0.0"));
	lines.push_back(string("ATOM       2      E1          EO        44.0526    0.00   0.0"));
	int size=lines.size();

	vector<atom> atoms;
	for (int i=0; i<size; i++)
	{
		atoms.push_back(atom());
		atoms[i].readAtomLine(lines[i]);
		cout<<setw(10)<<setprecision(10) << atoms[i].atomOutputIndex(5) ;
		cout<<setw(10)<<atoms[i].atomOutputName() ;
		cout<<setw(10)<<atoms[i].atomOutputType() ;
		cout<<setw(10)<<setprecision(4)<<setiosflags(ios::fixed)<<atoms[i].atomOutputMass() ;
		cout<<setw(10)<<setprecision(2)<<setiosflags(ios::fixed)<<atoms[i].atomOutputCharge() ;
		cout<<setw(10)<<setprecision(1)<<setiosflags(ios::fixed)<<atoms[i].atomOutputUnset()<<endl ;
		atoms[i].atomOutputAtom(5);
		
	}
}*/
