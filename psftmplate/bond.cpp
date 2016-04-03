#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

class bond 
{
public:
	int bondIndex;
	int bondType;
	int bondList[2];
	string bondAtomname[2];
	bond()
	{
	bondIndex=0;
	bondType=0;
	bondList[0]=0, bondList[1]=0;
	bondAtomname[0]="";
	bondAtomname[1]="";
	}
	void readBondLine(const string &text);
	int bondOutputIndex(int shift);
	int bondOutputType();
	int bondOutputList0(int atomshift0);
	int bondOutputList1(int atomshift1);
	void bondOutputList(int atomshift0, int atomshift1);
	const string bondOutputAtomname0();
	const string bondOutputAtomname1();
	void bondOutputAtomname();
	void  bondOutputBond(int indexshit, int atomshift0, int atomshift1);
	~bond(){};
private:
};

void bond::readBondLine(const string &text)
{
	string bondTitle;
	//sscanf(text, "%s%d%d%d%d%s%s", &bondTitle, &bondIndex, &bondType, &bondList[0], &bondList[1], &bondAtomname[0], &bondAtomname[1]);
	istringstream str(text);
	str >> bondTitle >> bondIndex >> bondType >> bondList[0] >> bondList[1] >> bondAtomname[0] >> bondAtomname[1];

}

int bond::bondOutputIndex(int indexshift)
{
	//printf("%10d\n", bondIndex+indexshift);
	return(bondIndex+indexshift);
	//cout<<atomIndex+shift<<endl;
}

int bond::bondOutputType()
{
	//printf("%10d\n", bondType);
	return(bondType);
}

int bond::bondOutputList0(int atomshift0)
{	
	return(bondList[0]+atomshift0);
}

int bond::bondOutputList1(int atomshift1)
{
        return(bondList[1]+atomshift1);
}

void bond::bondOutputList(int atomshift0, int atomshift1)
{
	//printf("%10d\t%10d\n", bondOutputList0(atomshift0), bondOutputList1(atomshift1));
	cout <<setw(10) <<bondOutputList0(atomshift0);
	cout <<setw(10) <<bondOutputList1(atomshift1);
}

const string bond::bondOutputAtomname0()
{
	return(bondAtomname[0]);
}

const string bond::bondOutputAtomname1()
{
        return(bondAtomname[1]);
}
void bond::bondOutputAtomname()
{
	//printf("%10s\t%10s\n", bondOutputAtomname0(), bondOutputAtomname1());
	cout <<setw(10) <<bondOutputAtomname0();
	cout <<setw(10) <<bondOutputAtomname1();
}

void bond::bondOutputBond(int indexshift, int atomshift0, int atomshift1)
{
	//printf("%10d\t%10d\t%10d\t%10d\t%10s\t%10s\n", bondOutputIndex(indexshift), bondOutputType(), 
	//	bondOutputList0(atomshift0), bondOutputList1(atomshift1), bondOutputAtomname0(), bondOutputAtomname1());
	cout <<setw(10) <<bondOutputIndex(indexshift);
	cout <<setw(10) <<bondOutputType();
	cout <<setw(10) <<bondOutputList0(atomshift0);
	cout <<setw(10) <<bondOutputList1(atomshift1);
	cout <<setw(10) <<bondOutputAtomname0();
	cout <<setw(10) <<bondOutputAtomname1()<<endl;
}
/*int main(void)
{
	vector<string> lines;
	lines.push_back(string("BOND    1    1    1    2   OA  EO"));
	lines.push_back(string("BOND    2    2    2    3   EO  EO"));
	int size=lines.size();

	vector<bond> bonds;
	for (int i=0; i<size; i++)
	{
		bonds.push_back(bond());
		bonds[i].readBondLine(lines[i]);
		cout <<setw(10) <<bonds[i].bondOutputIndex(5);
		cout <<setw(10) <<bonds[i].bondOutputType();
		cout <<setw(10) <<bonds[i].bondOutputList0(5);
		cout <<setw(10) <<bonds[i].bondOutputList1(5);
		cout <<setw(10) <<bonds[i].bondOutputAtomname0();
		cout <<setw(10) <<bonds[i].bondOutputAtomname0()<<endl;
		bonds[i].bondOutputBond(5, 5, 5);
	}	
}*/
