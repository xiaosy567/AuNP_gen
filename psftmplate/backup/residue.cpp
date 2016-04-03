#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include "atom.cpp"
#include "bond.cpp"
#include "angle.cpp"

using namespace std;

class residue 
{
public:
	string residName;
	int  residNumAtom, residNumBond, residNumAngle, residNumIpr; 
	vector<atom> atoms;
	vector<bond> bonds;
	vector<angle> angles;
	residue()
	{
	residName="";
	residNumAtom=0, residNumBond=0, residNumAngle=0, residNumIpr=0;
	}
	void readResidFile(const string &fn);
	//Get the Residue Info
	string residGetName();
	int residGetNumAtom();
	int residGetNumBond();
	int residGetNumAngle();
	int residGetNumIpr();
	
	//Get Nth Atom Info
	int residGetNAtomIndex(int n, int shift);
	string residGetNAtomName(int n);
	string residGetNAtomType(int n);
	double residGetNAtomMass(int n);
	double residGetNAtomCharge(int n);
	double residGetNAtomUnset(int n);
	//Get Nth Bond info
	int residGetNBondIndex(int n, int shift);
	int residGetNBondType(int n);
	vector<int> residGetNBondAtomIndex(int n, int atomshift0, int atomshift1);
	vector<string> residGetNBondAtomName(int n);
	//Get Nth Angle Info	
	int residGetNAngleIndex(int n, int shift) ;
	int residGetNAngleType(int n)	;
	vector<int> residGetNAngleAtomIndex(int n, int atomshift0, int atomshift1, int atomshift2);
	vector<string> residGetNAngleAtomName(int n);
	
	~residue(){};
private:
};

void residue::readResidFile(const string &fn)
{
	cout<<"==READING RESIDUE INF FROM "<<fn<<"=="<<endl;	
	fstream inpfile(fn.c_str());
	
	if(!inpfile)
	{
		cout <<"Error!"<<"Failed In Open "<<fn<<endl;
		abort();
	}

	string line;
	string residTitle;
	string numAtomTitle;
	string numBondTitle;
	string numAngleTitle;
	string numImprTitle;
	int nResidFileTitle=5;
	//5 Title lines; 
	int nnline=0; 
	for (int i=0; i<nResidFileTitle; i++)
	{
		nnline++;
		if(!getline(inpfile, line))
		{
		 	cout <<"Error!" <<"Failed in read " << nnline <<"th line..."<<endl;
		}
		istringstream str(line);
		switch(i)
		{
		case(0):
			str >> residTitle >> residName ; break ;
		case(1):
			str >> numAtomTitle >> residNumAtom ; break ;
		case(2):
			str >> numBondTitle >> residNumBond ; break ;
		case(3):
			str >> numAngleTitle >> residNumAngle ; break ; 
		case(4):
			str >> numImprTitle >> residNumIpr ; break ; 
		} 
	}
	if (nnline != nResidFileTitle)
	{
		 cout <<"Error!" <<"Not Correct in Read the Title..." << endl;
	}
	
	//ATOM	
	for (int i=0; i<residNumAtom; i++)
	{
		nnline++;
		if(!getline(inpfile, line))
                {
                        cout <<"Error!" <<"Failed in read " << nnline <<"th line..."<<endl;
                }
                istringstream str(line);	
		atoms.push_back(atom());
		atoms[i].readAtomLine(line);
		cout<<"ATOM CHECK ||"; 
		atoms[i].atomOutputAtom(0);
	}
	if (nnline != (residNumAtom+nResidFileTitle))
        {
                 cout <<"Error!" <<"Not Correct in Read the ATOM..." << endl;
        }
	
	//BOND
	for (int i=0; i<residNumBond; i++)
	{
		nnline++;
		if(!getline(inpfile, line))
		{
			cout <<"Error!" <<"Failed in read " << nnline <<"th line..."<<endl;
		}
		istringstream str(line);
		bonds.push_back(bond());
		bonds[i].readBondLine(line);
		cout <<"BOND CHECK ||";
		bonds[i].bondOutputBond(0, 0, 0) ;
		
	}
	if (nnline != (residNumBond + residNumAtom + nResidFileTitle))
	{
		cout <<"Error!" <<"Not Correct in Read the BOND..." << endl;
	}
	
	//ANGLE
	for (int i=0; i<residNumAngle; i++)
	{
		nnline++;
		if(!getline(inpfile, line) )
		{
			cout <<"Error!" <<"Failed in read " << nnline <<"th line..."<<endl;
		}
		istringstream str(line);
		angles.push_back(angle());
		angles[i].readAngleLine(line);
		cout<<"ANGLE CHECK||";
		angles[i].angleOutputAngle(0, 0, 0, 0) ;
	}
	if (nnline != (residNumAngle+ residNumBond + residNumAtom + nResidFileTitle))	
	{
		cout <<"Error!" <<"Not Correct in Read the ANGLE..." <<endl;
	}
	
}

//Get the Residue Info
string residue::residGetName()
{	
	return(residName);
}

int residue::residGetNumAtom()
{
	return(residNumAtom);	
}

int residue::residGetNumBond()
{
	return(residNumBond);
}

int residue::residGetNumAngle()
{
	return(residNumAngle);
}

int residue::residGetNumIpr()
{
	return(residNumIpr);
}

//Get Nth Atom Info
int residue::residGetNAtomIndex(int n, int shift)
{
	return(atoms[n].atomOutputIndex(shift));
}

string residue::residGetNAtomName(int n)
{
	return(atoms[n].atomOutputName());
}

string residue::residGetNAtomType(int n)
{
	return(atoms[n].atomOutputType());
}

double residue::residGetNAtomMass(int n)
{
	return(atoms[n].atomOutputMass());
}

double residue::residGetNAtomCharge(int n)
{
	return(atoms[n].atomOutputCharge());
}

double residue::residGetNAtomUnset(int n)
{
	return(atoms[n].atomOutputUnset());
}

//Get the Nth Bond Info
int residue::residGetNBondIndex(int n, int shift)
{
	return(bonds[n].bondOutputIndex(shift));
}

int residue::residGetNBondType(int n)
{
	return(bonds[n].bondOutputType());
}

vector<int> residue::residGetNBondAtomIndex(int n, int atomshift0, int atomshift1)
{
	vector<int> bondAtomIndex(2);
	bondAtomIndex[0]=bonds[n].bondOutputList0(atomshift0);
	bondAtomIndex[1]=bonds[n].bondOutputList1(atomshift1);
	return(bondAtomIndex);
}

vector<string> residue::residGetNBondAtomName(int n)
{
	vector<string> bondAtomName(2);
	bondAtomName[0]=bonds[n].bondOutputAtomname0();
	bondAtomName[1]=bonds[n].bondOutputAtomname1();
	return(bondAtomName);
}

//Get the Nth Angle Info
int residue::residGetNAngleIndex(int n, int shift)
{
	return(angles[n].angleOutputIndex(shift));
}

int residue::residGetNAngleType(int n)
{
	return(angles[n].angleOutputType());
}

vector<int> residue::residGetNAngleAtomIndex(int n, int atomshift0, int atomshift1, int atomshift2)
{
	vector<int> angleatomindex(3);
	angleatomindex[0]=angles[n].angleOutputList0(atomshift0) ;
	angleatomindex[1]=angles[n].angleOutputList1(atomshift1) ;
	angleatomindex[2]=angles[n].angleOutputList2(atomshift2) ;
	return(angleatomindex);
}

vector<string> residue::residGetNAngleAtomName(int n)
{
	vector<string> angleatomname(3);
	angleatomname[0]=angles[n].angleOutputAtomname0();
	angleatomname[1]=angles[n].angleOutputAtomname1();
	angleatomname[2]=angles[n].angleOutputAtomname2();
	return(angleatomname);
}

/*int main()
{
	residue aa;
	string file;
	for (int j=0; j<2; j++)
	{
	switch(j)
	{
	case(0): 
	{
	 	file="c6e6.dat"; break;
	}
	case(1):
		file="wat.dat"; break;	
	default:
		cout<< "ERROR:"<<endl; 
		abort();
	}
	aa.readResidFile(file)	;
	for (int i=0; i<aa.residGetNumAtom(); i++)
	{
		cout << setw(10)<<setprecision(10)<<aa.residGetNAtomIndex(i, 0);
		cout<<setw(10)<<aa.residGetNAtomName(i);
		cout<<setw(10)<<aa.residGetNAtomType(i);
		cout<<setw(10)<<setprecision(4)<<setiosflags(ios::fixed)<<aa.residGetNAtomMass(i);
		cout<<setw(10)<<setprecision(2)<<setiosflags(ios::fixed)<<aa.residGetNAtomCharge(i);
		cout<<setw(10)<<setprecision(1)<<setiosflags(ios::fixed)<<aa.residGetNAtomUnset(i)<<endl;;
	}
	vector< vector<int> > bondlist;
	vector< vector<string> > bondnamelist;
	for (int i=0; i<aa.residGetNumBond(); i++)
	{
		cout <<setw(10)<<aa.residGetNBondIndex(i, 0);
		cout <<setw(10)<<aa.residGetNBondType(i);
		bondlist.push_back(aa.residGetNBondAtomIndex(i, 0, 0));
		cout <<setw(10) <<bondlist[i][0];
		cout <<setw(10) <<bondlist[i][1];
		bondnamelist.push_back(aa.residGetNBondAtomName(i));
		cout <<setw(10) <<bondnamelist[i][0];
		cout <<setw(10) <<bondnamelist[i][1];
		cout<<"\n";
	}
	vector< vector <int> > anglelist;
	vector< vector <string> > anglenamelist;
	for (int i=0; i< aa.residGetNumAngle(); i++)
	{
		cout <<setw(10)<<aa.residGetNAngleIndex(i, 0);
		cout <<setw(10)<<aa.residGetNAngleType(i);
		anglelist.push_back(aa.residGetNAngleAtomIndex(i, 0, 0, 0));
		cout <<setw(10)<<anglelist[i][0];
		cout <<setw(10)<<anglelist[i][1];
		cout <<setw(10)<<anglelist[i][2];
		anglenamelist.push_back(aa.residGetNAngleAtomName(i));
		cout <<setw(10)<<anglenamelist[i][0];
		cout <<setw(10)<<anglenamelist[i][1];
		cout <<setw(10)<<anglenamelist[i][2];
		cout <<"\n";
	}
	}
}*/
