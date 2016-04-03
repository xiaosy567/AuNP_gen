#include <iostream>
#include <vector>
#include <utility>
using namespace std;

int main()
{
vector<vector<int> > bonds;
vector<int> bond(2);
bond[0]=1;
bond[1]=2;
bonds.push_back(bond);

for (int i=0; i<bonds.size(); i++)
{
 cout << bonds[i][0] << "-->" << bonds[i][1] << endl;
}

vector< vector<int> > angles;
vector<int> angle(3);
angle[0]=1;
angle[1]=2;
angle[2]=3;
angles.push_back(angle);

for (int i=0; i<angles.size(); i++)
{
 cout << angles[i][0] << "-->" 
	<< angles[i][1] << "-->"
	angles[i][2] << endl;
}

return 0;
}
