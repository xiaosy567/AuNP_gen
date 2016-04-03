#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "residue.cpp"
#include "function.h"

using namespace std;

class sys
{
public:
	vector<residue> molResids;            //each residue type will add to "molResidList" vector;
	vector<int> molResidAmount;		//the amount of each residue type;

	vector< vector<int> > molResidIndexIF;  //ResidIndexIF: save the the inital and final index number
						//molResidIndexIF[i][0]: the i residue's inital index number
						//molResidIndexIF[i][1]: the i residue's final index number
						
	vector<string> molResidPsfList;		//the name of each residue's PSF file will be add to ""
						//"molResidPsfList" vector
	vector<string> AtomTypeList;	//Atom type list;
	vector< vector<string> > BondTypeList ;	//Bond type list;	
	vector< vector<string> > AngleTypeList; //Angle Type list;
	
	sys() {};
	void readSysFile(const string &fn);

	//Generate/output the atom type list
	void sysGenAtomTypeList();
	void sysOutputAtomTypeList();
	//Generate/output the bond type list
	void sysGenBondTypeList();
	void sysOutputBondTypeList();
	//Generate/output the angle type list
	void sysGenAngleTypeList();
	void sysOutputAngleTypeList();

	//Write atom/bond/angle lists to file
	void sysWriteTypeListToFile();	

	//write psf 
	void sysWritePSF();
	//write lammps data file
	void sysWriteLammpsDataFile();

	//Write ALL
	void sysWriteAll();
	

	~sys(){};
private:
};

void sys::readSysFile(const string &fn)
{
	cout<<setw(25)<<""<<setw(28)<<"==READING SYSTEM INFO FROM  "<<fn<<"=="<<endl;
	fstream inpfile(fn.c_str());
	if(!inpfile)
	{
		cout<<"Error!"<<"Failed In Open "<<fn<<endl;
		abort();
	}
	string line;
	string residName;
	int residAmount;
	int residIndexi;
	int residIndexf;
	string residPsfFile;
        vector<int> ResidIndexIF(2);            //of each residue type; Then this vector will be add to 
        string firstchar;  
        string lastchar;   	                                //the "molResidIndexIF" vector;
        int lenstr=0;
        int i=0; 
	while(getline(inpfile, line))
	{
		//lenstr=line.strlen();
		lenstr=line.length();
                if( lenstr<2 ) { continue; }
                else {
			firstchar=line.substr(0,1);
			lastchar=line.substr(lenstr-1, 1);
			if( firstchar==string("#") ||
                    	lastchar!=string(";"))	
                	{
				continue;
			}
		}
		istringstream str(line);
		str >> residName >> residAmount >> residIndexi >> residIndexf >> residPsfFile;
		//If the last charactor of residPsfFile is ";", erase it;
		//cout<<residPsfFile<<endl;
		if(lastchar==string(";"))
		{
			residPsfFile.erase(residPsfFile.length()-1, 1);
		}
		//cout<<residPsfFile<<endl;

		//save the ith residue type
		molResids.push_back(residue());
		molResids[i].readResidFile(residPsfFile);
		//save the ith residue type's amount
		molResidAmount.push_back(residAmount);
		//save the ith residue type index (initial: final)
		ResidIndexIF[0]=residIndexi;
		ResidIndexIF[1]=residIndexf;
		//save the ith residue type's PsfFile name
		molResidPsfList.push_back(residPsfFile);
		i++;
	}

	//Generate atom type list
	sysGenAtomTypeList();
	//Generate bond type list
	sysGenBondTypeList();
	//Generate angle type list
	sysGenAngleTypeList();

	cout<<setw(25)<<""<<setw(28)<<"==END READ SYSTEM INFO FROM "<<fn<<"=="<<endl;
	
	inpfile.close();
}

//Generate/output atom type list
void sys::sysGenAtomTypeList()
{
	string atomName;
	string temp;
	bool   point;	
	for(int i=0; i< molResids.size(); i++)
	{
		for(int j=0; j< molResids[i].residGetNumAtom(); j++)
		{
			atomName=molResids[i].residGetNAtomType(j);
			if(i==0 && j==0)
			{
				AtomTypeList.push_back(atomName);
			}
			else
			{
				for(int k=0; k<AtomTypeList.size(); k++)
				{
					temp=AtomTypeList[k];
					point=check_strings(temp, atomName);
					if(point){continue;}
				}
				if(!point){AtomTypeList.push_back(atomName);}
			}
		}
	}
}

void sys::sysOutputAtomTypeList()
{
	for(int i=0; i<AtomTypeList.size(); i++)
	{
		cout<<i+1<<"  "<<AtomTypeList[i]<<endl;
	}
}

//Generate/output bond type list
void sys::sysGenBondTypeList()
{
	string bondatom1, bondatom2, bondatom3, bondatom4;
	vector<string> bondatom_pair(2);
	string temp;
	bool point=false;
	for(int i=0; i< molResids.size(); i++)
	{
		for(int j=0; j< molResids[i].residGetNumBond(); j++)
		{
			bondatom_pair=molResids[i].residGetNBondAtomName(j);
			
			//cout<<bondatom_pair[0] << "  " <<bondatom_pair[1]<<endl;

			if(i==0 && j==0)
			{
				BondTypeList.push_back(bondatom_pair);
			}
			else
			{
					bondatom1=bondatom_pair[0];
					bondatom2=bondatom_pair[1];
				for(int k=0; k<BondTypeList.size(); k++)
				{
					bondatom3=BondTypeList[k][0];
					bondatom4=BondTypeList[k][1];
					
					point=check_strings(bondatom1, bondatom2, bondatom3, bondatom4);
					if(point){break;}
				}
				if(!point)
				{
					BondTypeList.push_back(bondatom_pair);
				}
			}
		}
	}
}


void sys::sysOutputBondTypeList()
{
	for(int i=0; i<BondTypeList.size(); i++)
	{
		 cout<<i+1<<"  "<<BondTypeList[i][0]<<"  "<<BondTypeList[i][1]<<endl;
	}
	
}

//Generate/output Angle type list
void sys::sysGenAngleTypeList()
{
        string angleatom1, angleatom2, angleatom3, angleatom4, angleatom5, angleatom6;
        vector< vector<string> > angleatoms;
        vector<string> angleatom_triple(3);
        string temp;
        bool point;
        for(int i=0; i< molResids.size(); i++)
        {   
                for(int j=0; j< molResids[i].residGetNumAngle(); j++)
                {   
			//cout<<"debug:"<<"j="<<j<<endl;
                        //angleatoms.push_back(molResids[i].residGetNAngleAtomName(j));
                        //angleatom_triple[0]=angleatoms[j][0];
                        //angleatom_triple[1]=angleatoms[j][1];
			//angleatom_triple[2]=angleatoms[j][2];
			angleatom_triple=molResids[i].residGetNAngleAtomName(j);

                        if(i==0 && j==0)
                        {   
                                AngleTypeList.push_back(angleatom_triple);
				//sysOutputAngleTypeList();
                        }   
                        else
                        {   
					point=true;
                                        //angleatom1=angleatoms[j][0];
                                        //angleatom2=angleatoms[j][1];
					//angleatom3=angleatoms[j][2];
					angleatom1=angleatom_triple[0];
					angleatom2=angleatom_triple[1];
					angleatom3=angleatom_triple[2];
                                for(int k=0; k<AngleTypeList.size(); k++)
                                {   
                                        angleatom4=AngleTypeList[k][0];
                                        angleatom5=AngleTypeList[k][1];
					angleatom6=AngleTypeList[k][2];
				
					//cout<<"debug:"<<"k="<<k<<endl;
    
                                        point=check_strings(angleatom1, angleatom2, angleatom3, angleatom4, angleatom5, angleatom6);
					//cout<<"debug:"<<"point="<<point<<endl;
                                        if(point){break;}
					//cout<<"debug:"<<"continue"<<endl;
                                }   
				//cout<<"before ouput!"<<endl;
                                if(!point)
                                {   
					//cout<<"debug:"<<"AngleTypeList.push_back(angleatom_triple)"<<endl;
                                        AngleTypeList.push_back(angleatom_triple);
					//cout<<"debug:"<<"AngleTypeList.push_back(angleatom_triple)..."<<endl;
					//sysOutputAngleTypeList();
                                }   
                        }   
                }   
        }   
}

void sys::sysOutputAngleTypeList()
{
        for(int i=0; i<AngleTypeList.size(); i++)
        {
                 cout<<i+1<<"  "<<AngleTypeList[i][0]<<"  "<<AngleTypeList[i][1]<<"  "<<AngleTypeList[i][2]<<endl;
        }

}

void sys::sysWritePSF()
{
	cout<<"==Write the PSF file of the defined system=="<<endl;
	int NATOM=0, NBOND=0, NTHETA=0;
	for (int i=0; i<molResids.size(); i++)
	{
	NATOM = NATOM  + molResids[i].residGetNumAtom()*molResidAmount[i];
	NBOND = NBOND  + molResids[i].residGetNumBond()*molResidAmount[i];
	NTHETA= NTHETA + molResids[i].residGetNumAngle()*molResidAmount[i];
	}

	string psfFileName="OUT.psf";
	fstream outfile(psfFileName.c_str(), ios::out);
	if(!outfile)
	{
		cout<<"==Failed in OPEN/CREATE file: OUT.psf=="<<endl;
	}
	//PSF file title
	cout<<"PSF"<<endl;	
	outfile<<"PSF"<<endl;
	//cout<<""<<endl;
	outfile<<""<<endl;
	//cout<<setw(7)<<"2"<<setw(1)<<" "<<"!NTITLE"<<endl;
	outfile<<setw(7)<<"2"<<" "<<"!NTITLE"<<endl;
	//cout<<"REMARKS"<<endl;
	outfile<<"REMARKS"<<endl;
	//cout<<"REMARKS"<<endl;
	outfile<<"REMARKS"<<endl;
	//cout<<""<<endl;
	outfile<<""<<endl;
	//
	//Write ATOMS
	cout<<setw(8)<<NATOM<<" "<<"!NATOM"<<endl;
	outfile<<setw(8)<<NATOM<<" "<<"!NATOM"<<endl;
	int i, j, k;
	int NatomPerResidue=0;
	int AtomIndex=0;
	int ResidueIndex=0;
	string ResidueName="";
	for (i=0; i<molResids.size(); i++)
	{
		NatomPerResidue = molResids[i].residGetNumAtom();	
		ResidueName = molResids[i].residGetName();
		for(j=0; j<molResidAmount[i]; j++)
		{
			ResidueIndex = ResidueIndex + 1;
			for(k=0; k<molResids[i].residGetNumAtom(); k++)
			{
				AtomIndex = AtomIndex + 1;
				outfile<<setw(8)<<AtomIndex;
				outfile<<setw(1)<<" ";
				outfile<<setw(4)<<ResidueName;
				outfile<<setw(5)<<ResidueIndex;
				outfile<<setw(1)<<" ";
				outfile<<setw(4)<<ResidueName;
				outfile<<setw(1)<<" ";
				outfile<<setw(4)<<molResids[i].residGetNAtomName(k);
				outfile<<setw(1)<<" ";
				outfile<<setw(4)<<molResids[i].residGetNAtomType(k);
				outfile<<setw(1)<<"";
				outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<molResids[i].residGetNAtomCharge(k);
				outfile<<setw(7)<<"";
				outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<molResids[i].residGetNAtomMass(k);
				outfile<<"\n";
			}
		}
	}
	if (AtomIndex==NATOM)
	{
		cout<<"==The ATOM PART OF PSF FILE WAS WRITTEN CORRECTLY=="<<endl;
		outfile<<""<<endl;
	}
	
	//Write BOND
	cout<<setw(8)<<NBOND<<setw(1)<<" "<<"!NBOND"<<endl;
	outfile<<setw(8)<<NBOND<<setw(1)<<" "<<"!NBOND"<<endl;
	int AtomIndexShift=0; 
	int linecontrol=0;
	vector<int> bondatom_pair(2);
	for (i=0; i<molResids.size(); i++)
	{
		for(j=0; j<molResidAmount[i]; j++)
		{
			for(k=0; k<molResids[i].residGetNumBond(); k++)
			{
				linecontrol = linecontrol + 1;
				bondatom_pair = molResids[i].residGetNBondAtomIndex(k, AtomIndexShift, AtomIndexShift);
				outfile<<setw(8)<<bondatom_pair[0]<<setw(8)<<bondatom_pair[1];
				if((linecontrol%4)==0) 
				{
					outfile<<"\n";
				}
				/*if(k==(molResids[i].residGetNumBond()-1)) 
				{
					AtomIndexShift = AtomIndexShift + molResids[i].residGetNumAtom();
				}*/
			}
			AtomIndexShift = AtomIndexShift + molResids[i].residGetNumAtom();
		}
	}
	if (linecontrol==NBOND) 
	{
		cout<<"==The BOND PART OF PSF FILE WAS WRITTEN CORRECTLY=="<<endl;
		outfile<<"\n"<<endl;
	}
	//Write ANGLE
	cout<<setw(8)<<NTHETA<<setw(1)<<" "<<"!NANGLE"<<endl;
	outfile<<setw(8)<<NTHETA<<setw(1)<<" "<<"!NANGLE"<<endl;
	AtomIndexShift=0;
	linecontrol=0;
	vector<int> angleatom_triple(3);
	for (i=0; i<molResids.size(); i++)
	{
		for(j=0; j<molResidAmount[i]; j++)
		{
			for(k=0; k<molResids[i].residGetNumAngle(); k++)
			{
				linecontrol = linecontrol + 1;
				angleatom_triple=molResids[i].residGetNAngleAtomIndex(k, AtomIndexShift, AtomIndexShift, AtomIndexShift);
				outfile<<setw(8)<<angleatom_triple[0]<<setw(8)<<angleatom_triple[1]<<setw(8)<<angleatom_triple[2];
				if((linecontrol%3)==0)
				{
					outfile<<"\n";
				}
			}
			AtomIndexShift = AtomIndexShift + molResids[i].residGetNumAtom();
		}
	}	
	if (linecontrol == NTHETA )
	{
		cout<<"==The ANGLE PART OF PSF FILE WAS WRITTEN CORRECTLY=="<<endl;
	}
}

//Write the lammmps formate data file, however, no coordinates is included;
void sys::sysWriteLammpsDataFile()
{
	cout<<"==Write the Lammps Data File of the defined system=="<<endl;
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	cout<<setw(5)<<""<<"|          LAMMPS DATA FILE        |"<<"\n";
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	int NATOM=0, NBOND=0, NTHETA=0;
	for (int i=0; i<molResids.size(); i++)
	{
		NATOM = NATOM  + molResids[i].residGetNumAtom()*molResidAmount[i];
		NBOND = NBOND  + molResids[i].residGetNumBond()*molResidAmount[i];
		NTHETA= NTHETA + molResids[i].residGetNumAngle()*molResidAmount[i];
	}
	
	string lammpsDataFilenName="OUT.lammps.dat";
	fstream outfile(lammpsDataFilenName.c_str(), ios::out);
	if(!outfile)
	{
		cout<<"==Failed in OPEN/CREATE file: "<<lammpsDataFilenName<<"=="<<endl;
	}
	
	//lammps data file title;
	cout<<"LAMMPS Description"<<endl;
	outfile<<"LAMMPS Description"<<endl;
	cout<<"\n";
	outfile<<"\n";
	//set the amount of atom/bond/angle;
	outfile<<setw(12)<<NATOM<<setw(2)<<""<<"atoms"<<"\n";
	outfile<<setw(12)<<NBOND<<setw(2)<<""<<"bonds"<<"\n";
	outfile<<setw(12)<<NTHETA<<setw(2)<<""<<"angles"<<"\n";
	outfile<<"\n";
	//set the number of atom/bond/angle type
	outfile<<setw(12)<<AtomTypeList.size()<<setw(2)<<""<<"atom types"<<"\n";
	outfile<<setw(12)<<BondTypeList.size()<<setw(2)<<""<<"bond types"<<"\n";
	outfile<<setw(12)<<AngleTypeList.size()<<setw(2)<<""<<"angle types"<<"\n";
	outfile<<"\n";
	//set the box
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(11)<<"xlo xhi"<<"\n";
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(11)<<"ylo yhi"<<"\n";
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(9)<<setprecision(4)<<setiosflags(ios::fixed)<<0.0;
	outfile<<setw(11)<<"zlo zhi"<<"\n";
	outfile<<"\n";

	//ATOMS PART
	cout<<"==Write the ATOM part of lammps data file=="<<endl;
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	cout<<setw(5)<<""<<"|Atoms line syntax: molecular/angle|"<<"\n";
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	outfile<<"Atoms"<<endl;
	outfile<<"\n";
	
	int i, j, k, m;
	int NatomPerResidue=0;
	int AtomIndex=0;
	int ResidueIndex=0;
	string AtomName="", TypeName="";
	int AtomTypeIndex=0;
	bool typecontrol=false;
	//Charge bool
	bool boolCharge=false;
	
	//If the atoms are all neutral, the "atom style" 
	//will be "atomic", otherwise, the outut style will be "charge"
	//
	for (i=0; i<molResids.size(); i++)
	{
		for(j=0; j<molResidAmount[i]; j++)
		{
			for(k=0; k<molResids[i].residGetNumAtom(); k++)
			{
				if(fabs(molResids[i].residGetNAtomCharge(k)-0.0) <= 0.001 )
				{
					boolCharge=true;
					goto JUMP;
				}
			}
		}
	}	
	
 	JUMP:;

	for (i=0; i<molResids.size(); i++)
	{
		NatomPerResidue = molResids[i].residGetNumAtom();
		for(j=0; j<molResidAmount[i]; j++)
		{
			ResidueIndex = ResidueIndex + 1;
			for(k=0; k<molResids[i].residGetNumAtom(); k++)
			{
				AtomIndex = AtomIndex + 1;
				AtomName = molResids[i].residGetNAtomType(k);
				typecontrol = false;
				//cout<<"i="<<i<<"  "<<"j="<<j<<"  "<<"k="<<k<<endl;
				for (m=0; m<AtomTypeList.size(); m++)
				{
					TypeName=AtomTypeList[m];
					typecontrol=check_strings(TypeName, AtomName);
					//cout<<AtomName<<"  "<<TypeName<<endl;
					if(typecontrol)
					{
						AtomTypeIndex = m + 1;
						//cout<<"in m circle:typecontrol="<<typecontrol<<endl;
						break;
					}	
				}
				if(!typecontrol)
				{
					cout<<"=="<<"The "<<AtomName<<" was not found in the AtomgTypeList=="<<endl;
					//cout<<"typecontrol="<<typecontrol<<endl;
				}
				outfile<<setw(8)<<AtomIndex;
				outfile<<setw(1)<<" ";
				outfile<<setw(5)<<ResidueIndex;
				outfile<<setw(1)<<"";
				outfile<<setw(4)<<AtomTypeIndex;
				if(boolCharge)
				{
					outfile<<setw(9)<<setprecision(5)<<setiosflags(ios::fixed)<<molResids[i].residGetNAtomCharge(k);
				}
				outfile<<setw(9)<<setprecision(3)<<setiosflags(ios::fixed)<<0.0;
				outfile<<setw(1)<<"";
				outfile<<setw(9)<<setprecision(3)<<setiosflags(ios::fixed)<<0.0;
				outfile<<setw(1)<<"";
				outfile<<setw(9)<<setprecision(3)<<setiosflags(ios::fixed)<<0.0;
				outfile<<"\n";
			}
		}
	}
	if (AtomIndex==NATOM)
	{
		cout<<"==Write the ATOM part of lammps data file CORRECTLY=="<<endl;
		outfile<<"\n";
	}
	
	//BOND
	cout<<"==Write the BOND part of lammps data file=="<<endl;
	outfile<<"\n";
	outfile<<"Bonds"<<"\n";
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	cout<<setw(5)<<""<<"|              Bonds               |"<<"\n";
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	outfile<<"\n";
	int AtomIndexShift=0;
	int bond_index=0;
	int BondTypeIndex=0;
	vector<int> bondatom_pair_index(2);
	vector<string> bondatom_pair_name(2);
	string name1="", name2="", TypeName1="", TypeName2="";
	for (i=0; i<molResids.size(); i++)
	{
		for(j=0; j<molResidAmount[i]; j++)
		{
			for(k=0; k<molResids[i].residGetNumBond(); k++)
			{
				bond_index = bond_index + 1;
				bondatom_pair_index = molResids[i].residGetNBondAtomIndex(k, AtomIndexShift, AtomIndexShift);
				bondatom_pair_name = molResids[i].residGetNBondAtomName(k);
				name1 = bondatom_pair_name[0];
				name2 = bondatom_pair_name[1];
				typecontrol = false;
				BondTypeIndex=0;
				for (m=0; m<BondTypeList.size(); m++)
				{
						TypeName1 = BondTypeList[m][0];
						TypeName2 = BondTypeList[m][1];
						typecontrol = check_strings(bondatom_pair_name[0], bondatom_pair_name[1],
								BondTypeList[m][0], BondTypeList[m][1]);
						//typecontrol = check_strings(name1, name2, TypeName1, TypeName2);

						if(typecontrol)
						{
							BondTypeIndex =  m + 1;
							break;	
						}
				}
				outfile<<setw(8)<<bond_index;
				outfile<<setw(1)<<"";
				outfile<<setw(5)<<BondTypeIndex;
				outfile<<setw(1)<<"";
				outfile<<setw(8)<<bondatom_pair_index[0]<<setw(1)<<""<<setw(8)<<bondatom_pair_index[1]<<endl;
			}
			AtomIndexShift = AtomIndexShift + molResids[i].residGetNumAtom();
		}
	}
	
	if((AtomIndexShift==NATOM) && (bond_index == NBOND))
	{
		cout<<"==Write the BOND part of lammps data file CORRECTLY=="<<endl;
		outfile<<"\n";
	}

	//ANGLE
	cout<<"==Write the ANGLE part of lammps data file=="<<endl;
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";
	cout<<setw(5)<<""<<"|               ANGLE              |"<<"\n";
	cout<<setw(5)<<""<<"------------------------------------"<<"\n";

	outfile<<"Angles"<<"\n";
	outfile<<"\n";
	
	int angle_index=0;
	int AngleTypeIndex=0;
	vector<string> angleatom_triple_name(3);
	vector<int> angleatom_triple_index(3);	
	string name3="";
	AtomIndexShift=0;
	for (i=0; i<molResids.size(); i++)
	{
		for(j=0; j<molResidAmount[i]; j++)
		{
			for(k=0; k<molResids[i].residGetNumAngle(); k++)
			{
					angle_index = angle_index + 1;
					angleatom_triple_index = molResids[i].residGetNAngleAtomIndex(k, 
							AtomIndexShift, AtomIndexShift, AtomIndexShift);
					angleatom_triple_name = molResids[i].residGetNAngleAtomName(k);
	
					typecontrol = false;
					
					for (m=0; m<AngleTypeList.size(); m++)
					{
						typecontrol = check_strings(angleatom_triple_name[0], angleatom_triple_name[1],
								angleatom_triple_name[2], AngleTypeList[m][0],
								AngleTypeList[m][1], AngleTypeList[m][2]);
					
						if(typecontrol)
						{
							AngleTypeIndex = m + 1;
							break;
						}
					}
					outfile<<setw(8)<<angle_index;
					outfile<<setw(1)<<"";
					outfile<<setw(5)<<AngleTypeIndex;
					outfile<<setw(8)<<angleatom_triple_index[0]<<setw(8)<<angleatom_triple_index[1]
								<<setw(8)<<angleatom_triple_index[2]<<endl;
			}
			AtomIndexShift = AtomIndexShift + molResids[i].residGetNumAtom();
		}
	}
	
	if((AtomIndexShift==NATOM) && (angle_index == NTHETA))
	{
		cout<<"==Write the ANGLE part of lammps data file CORRECTLY=="<<endl;	
		outfile<<"\n";
	}
	

	//end
	cout<<"==Write Lammps Data File of defined system SUCCESSFULLY=="<<endl;
}

void sys::sysWriteTypeListToFile()
{
	cout<<"==Write the ATOM/BOND/ANGLE LISTS of the defined system=="<<endl;	
	
	string listsDataFilenName="OUT.list.lammps.dat";
	fstream outfile(listsDataFilenName.c_str(), ios::out);
	if(!outfile)
	{
		cout<<"==Failed in OPEN/CREATE file: "<<listsDataFilenName<<"=="<<endl;
	}
	
	//lists file title;
	outfile<<"#ATOM/BOND/ANGLE lists used for generating force field;"<<endl;
	outfile<<"\n";
	outfile<<"pair_style      cg/cmm 15.0"<<endl;
	outfile<<"bond_style      harmonic"<<endl;
	outfile<<"angle_style     harmonic"<<endl;
	outfile<<"special_bonds  lj       0.0 0.0 1.0"<<endl;;
	outfile<<"\n";

	//ATOM LIST
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";
        cout<<setw(5)<<""<<"|               ATOM               |"<<"\n";
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";

	int i=0;
	for (i=0; i<AtomTypeList.size(); i++)
	{
		outfile<<setw(4)<<"mass"<<setw(5)<<i+1<<setw(2)<<""<<setw(5)<<""<<setw(1)<<"#"<<setw(5)<<AtomTypeList[i]<<endl;
	}
	outfile<<"\n";
	
	//BOND LIST
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";
        cout<<setw(5)<<""<<"|              Bonds               |"<<"\n";
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";

	for (i=0; i<BondTypeList.size(); i++)
	{
		outfile<<"bond_coeff"<<setw(5)<<i+1<<setw(5)<<"";
		outfile<<setw(1)<<"#"<<setw(5)<<BondTypeList[i][0]<<setw(5)<<BondTypeList[i][1]<<endl;
	}
	outfile<<"\n";
	

	//ANGLE LIST
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";
        cout<<setw(5)<<""<<"|               ANGLE              |"<<"\n";
        cout<<setw(5)<<""<<"------------------------------------"<<"\n";

	for (i=0; i<AngleTypeList.size(); i++)
	{
		outfile<<"angle_coeff"<<setw(5)<<i+1<<setw(5)<<"";
		outfile<<setw(1)<<"#";
		outfile<<setw(5)<<AngleTypeList[i][0]<<setw(5)<<AngleTypeList[i][1]<<setw(5)<<AngleTypeList[i][2]<<endl;
	}

	//END
	cout<<"==Write ATOM/BOND/ANGLE LISTS  of  system  SUCCESSFULLY=="<<endl;
}

void sys::sysWriteAll()
{
	sysWritePSF();
	sysWriteLammpsDataFile();
	sysWriteTypeListToFile();
}
int main()
{	
	/*string line5="  USTC  ";
	//cout<<removeBEblanket(line5)<<endl;
	string line6="  USTC   ";
	//cout<<removeBEblanket(line6)<<endl;
	cout<<check_strings(line5, line6)<<endl;

	string line1="OA";
	string line2="OA";
	if (line1==line2)
	{
		cout<<"line1==line2"<<endl;
	}

	string line3="OA";
	string line4="OE";
	if ((line1==line3 && line2==line4) || (line1==line4 && line2==line3) )
	{
		cout<<"bond1==bond2"<<endl;
	}*/
	sys cc;
	string inputfile ;
	inputfile="INPUT" ;
	cc.readSysFile(inputfile) ;
	cc.sysOutputAtomTypeList();
	cc.sysOutputBondTypeList();
	cc.sysOutputAngleTypeList();
	
	cc.sysWriteAll();
	
	
}
