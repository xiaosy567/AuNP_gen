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

class angle 
{
public:
	int angleIndex;
	int angleType;
	int angleList[3];
	string angleAtomname[3];
	angle()
	{
	angleIndex=0;
	angleType=0;
	angleList[0]=0, angleList[1]=0, angleList[2]=0;
	angleAtomname[0]="";
	angleAtomname[1]="";
	angleAtomname[2]="";
	}
	void readAngleLine(const string &text);
	int angleOutputIndex(int indexshift);
	int angleOutputType();
	int angleOutputList0(int atomshift0);
	int angleOutputList1(int atomshift1);
	int angleOutputList2(int atomshift2);
	void angleOutputList(int atomshift0, int atomshift1, int atomshift2);
	const string angleOutputAtomname0();
	const string angleOutputAtomname1();
	const string angleOutputAtomname2();
	void angleOutputAtomname();
	void angleOutputAngle(int indexshift, int atomshift0, int atomshift1, int atomshift2);
	~angle(){};
private:
};

void angle::readAngleLine(const string &text)
{
	string angleTitle;
	//sscanf(text, "%s%d%d%d%d%d%s%s%s", &angleTitle, &angleIndex, &angleType, 
	//		&angleList[0], &angleList[1], &angleList[2], &angleAtomname[0], &angleAtomname[1], &angleAtomname[2]);
	istringstream str(text);
	str >> angleTitle >> angleIndex >> angleType >> angleList[0] >> angleList[1] >> angleList[2] >> angleAtomname[0] >> angleAtomname[1] >> angleAtomname[2];

}

int angle::angleOutputIndex(int indexshift)
{
	//printf("%10d\n", angleIndex+indexshift);
	return(angleIndex+indexshift);
	//cout<<atomIndex+shift<<endl;
}

int angle::angleOutputType()
{
	//printf("%10d\n", angleType);
	return(angleType);
}

int angle::angleOutputList0(int atomshift0)
{
	//printf("%10d", angleList[0]+atomshift0);
	return(angleList[0]+atomshift0);
}

int angle::angleOutputList1(int atomshift1)
{
        //printf("%10d", angleList[1]+atomshift1);
        return(angleList[1]+atomshift1);
}

int angle::angleOutputList2(int atomshift2)
{
        //printf("%10d", angleList[2]+atomshift2);
        return(angleList[2]+atomshift2);
}

void angle::angleOutputList(int atomshift0, int atomshift1, int atomshift2)
{
	//printf("%10d\t%10d\t%10d\n", angleOutputList0(atomshift0), angleOutputList1(atomshift1), angleOutputList2(atomshift2));
	cout << setw(10) << angleOutputList0(atomshift0);
	cout << setw(10) << angleOutputList1(atomshift1);
	cout << setw(10) << angleOutputList2(atomshift2);
}

const string angle::angleOutputAtomname0()
{
	return(angleAtomname[0]);
}

const string angle::angleOutputAtomname1()
{
	return(angleAtomname[1]);
}

const string angle::angleOutputAtomname2()
{
	return(angleAtomname[2]);
}

void angle::angleOutputAtomname()
{
	//printf("%10s\t%10s\t%10s\n", angleOutputAtomname0(), angleOutputAtomname1(), angleOutputAtomname2());
	cout <<setw(10) <<angleOutputAtomname0();
	cout <<setw(10) <<angleOutputAtomname1();
	cout <<setw(10) <<angleOutputAtomname2();
}

void angle::angleOutputAngle(int indexshift, int atomshift0, int atomshift1, int atomshift2)
{
	//printf("%10d\t%10d\t%10d\t%10d\t%10d\t%10s\t%10s\t%10s\n", angleOutputIndex(indexshift), angleOutputType(), 
	//	angleOutputList0(atomshift0), angleOutputList1(atomshift1), angleOutputList2(atomshift2),
	//	angleOutputAtomname0(), angleOutputAtomname1(), angleOutputAtomname2());
	cout << setw(10) << angleOutputIndex(indexshift);
	cout << setw(10) << angleOutputType();
	cout << setw(10) << angleOutputList0(atomshift0);
	cout << setw(10) << angleOutputList1(atomshift1);
	cout << setw(10) << angleOutputList2(atomshift2); 
	cout << setw(10) << angleOutputAtomname0();
	cout << setw(10) << angleOutputAtomname1(); 
	cout << setw(10) << angleOutputAtomname2()<<endl;
}

/*int main(void)
{
	vector<string> lines;
	lines.push_back(string("ANGLE   1    1    1    2   3   OA  EO  EO"));
	lines.push_back(string("ANGLE   2    2    2    3   4   EO  EO  EO"));
	int size=lines.size();
	
	vector<angle> angles;
	for (int i=0; i<size; i++)
	{
		angles.push_back(angle());
		angles[i].readAngleLine(lines[i]);
		cout << setw(10) << angles[i].angleOutputIndex(5);
		cout << setw(10) << angles[i].angleOutputType();
		cout << setw(10) << angles[i].angleOutputList0(5);
		cout << setw(10) << angles[i].angleOutputList1(5);
		cout << setw(10) << angles[i].angleOutputList2(5);
		cout << setw(10) << angles[i].angleOutputAtomname0();
		cout << setw(10) << angles[i].angleOutputAtomname1();
		cout << setw(10) << angles[i].angleOutputAtomname2()<<endl;
		angles[i].angleOutputAngle(5, 5, 5, 5);
	}
}*/
