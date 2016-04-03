#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

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
	atomMass=0.0;
	atomCharge=0.0;
	atomUnset=0.0;
	}
	void readline(const string &text);
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

void atom::readline(const string &text)
{
	string atomTitle;
	//sscanf(text, "%s%d%s%s%lf%lf%lf", atomTitle, &atomIndex, atomName, atomType, &atomMass, &atomCharge, &atomUnset);
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
	cout << atomOutputIndex(shift) << " "
		<< atomOutputName()<< " "
		<< atomOutputType()<< " "
		<< atomOutputMass()<< " "
		<< atomOutputCharge()<< " "
		<< atomOutputUnset()<<endl;	
}

int main(void)
{
	vector<string> lines;
	lines.push_back(string("ATOM       1      OA          OA        31.0337    0.00   0.0"));
	lines.push_back(string("ATOM       2      Ob          Oc        1.0337    0.00   0.0"));
	cout << lines.size() << endl;
	vector<atom> atoms;

	for(int i=0; i<2; i++)
	{
	atoms.push_back(atom());
	atoms[i].readline(lines[i]);
	cout << atoms[i].atomOutputIndex(5) << endl;
 	cout << atoms[i].atomOutputName() << endl;
	cout << atoms[i].atomOutputType() << endl;
	cout << atoms[i].atomOutputMass() << endl;
	cout << atoms[i].atomOutputCharge() << endl;
	cout << atoms[i].atomOutputUnset() << endl;
	atoms[i].atomOutputAtom(5);
	}
}
