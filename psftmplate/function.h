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

//remove the blanket the beginning and end of a string;
string removeBEblanket(const string &s);

//check the strings is the same or not;
bool check_strings(const string &s1, string &s2);  //s1 ?= s2
bool check_strings(const string &s1, string &s2, string &s3, string &s4); //(s1 s2) ?= (s3 s4)
bool check_strings(const string &s1, string &s2, string &s3, string &s4, string &s5, string &s6); // (s1 s2 s3) ?= (s4 s5 s6)

string removeBEblanket(const string &s) 
{
        string scopy=s;
        int scopySize=scopy.size();
	//cout<<"|"<<scopy<<"|"<<endl;
	/*for(int i=0; i<scopySize; i++)
	{
		cout<<"|"<<scopy.substr(i, 1)<<"|"<<endl;
	}*/
	//cout<<"- erase -"<<endl;
        while(scopy.substr(0, 1)==" ")
        {   
		//cout<<"|"<<scopy.substr(0, 1)<<"|"<<endl;
                scopy.erase(0, 1); //remove the first " ";
        }   
        while(scopy.substr(scopy.size()-1, 1)==" ")
        {   
                scopy.erase(scopy.size()-1, 1); // remove the last " ";
        }   
	//cout<<"|"<<scopy<<"|"<<endl;
        return(scopy);
}

//"check_strings" is used to check the three strings;
//They are same or not;
//such as: "OA"=="OA"
//         ("OA" "EO")==("EO" "OA") 
//         ("OA" "EO" "CT")==("CT" "EO" "OA")
bool check_strings(const string &s1, string &s2)
{
        string scopy1, scopy2;
        scopy1=removeBEblanket(s1);
        scopy2=removeBEblanket(s2);
        return(scopy1==scopy2);
}

bool check_strings(const string &s1, string &s2, string &s3, string &s4)
{
        string scopy1, scopy2, scopy3, scopy4;
        scopy1=removeBEblanket(s1);
        scopy2=removeBEblanket(s2);
        scopy3=removeBEblanket(s3);
        scopy4=removeBEblanket(s4);
        return((scopy1==scopy3 && scopy2==scopy4)||(scopy1==scopy4 && scopy2==scopy3) );
}


bool check_strings(const string &s1, string &s2, string &s3, string &s4, string &s5, string &s6)
{
        string scopy1, scopy2, scopy3, scopy4, scopy5, scopy6;
        scopy1=removeBEblanket(s1);
        scopy2=removeBEblanket(s2);
        scopy3=removeBEblanket(s3);
        scopy4=removeBEblanket(s4);
        scopy5=removeBEblanket(s5);
        scopy6=removeBEblanket(s6);
        return((scopy2==scopy5) &&( (scopy1==scopy4 && scopy3==scopy6)||(scopy1==scopy6 && scopy3==scopy4)) );
}

