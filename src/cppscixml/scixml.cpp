/** 
 * \file scixml.cpp
 * \brief Implements the scixml class
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: scixml.cpp 1454 2012-03-23 20:35:40Z gronoff $
 *
 */



#include "scixml.hpp"
using namespace std;

namespace fs = boost::filesystem;
#include <boost/version.hpp>

namespace ggfileutils
{

	std::string JoinPath(std::string vPath,std::string vFilename)
	{
		fs::path chem(vPath);
		fs::path file(vFilename);
		fs::path resu = chem/file;
		return resu.string();
	}



	std::string FilenameToPath(std::string vPathfile)
	{
		fs::path chem(vPathfile);
#if BOOST_VERSION <= 103500

		return chem.branch_path().string();
#else
		return chem.parent_path().string();
#endif
	}

	bool FileExists(std::string vPathfile)
	{
		fs::path chem(vPathfile);
		return fs::exists(chem);
	}

	std::string FullPath(std::string vPathfile)
	{
		fs::path chem(vPathfile);

		if(chem.has_root_directory())
			return vPathfile;
		fs::path pwd=fs::current_path();
		pwd/=chem;
		return pwd.string();
		//return chem.root_path().string();
	}

};
XmlParameters::XmlParameters(std::string vFilename)
{
	mFilename=vFilename;
	mMonteCarloActivated=false;
	mpDoc=new TiXmlDocument;
	if(!mpDoc->LoadFile(mFilename))
	{
		// It is put at debug because a real error is raised
		//Log::SetPriority(Log::DEBUGG,"Error raised  in XmlParameters");
		Log::mE<<"Impossible to load the xml document : "<<mFilename<<endl;
		Log::mE<<"Check if it is a valid xml document, check the path"<<endl;
		delete mpDoc;

		Error err("XmlParameters","Not valid xml document","Error while loading you document "+vFilename+". Please check if it exists and if it is xml valid (for example by opening it with firefox/iceweasel)");
		throw err;
		//exit(0);
	}
#ifdef DEBUG
	Log::mD<<"Votre document"<<mFilename<<" a ete charge"<<endl;
#endif
}


XmlParameters::~XmlParameters()
{
	delete mpDoc;
#ifdef DEBUG
	Log::mD<<"Fin de lecture de "<<mFilename<<", et destruction de l'objet associe"<<endl;
#endif
}



std::string XmlParameters::SubFileName(std::string vFilename)
{

	return JoinPath(trim(FilenameToPath(mFilename)),trim(vFilename));
}


std::string XmlParameters::GetSubFileName(std::string xPath)
{
	return SubFileName(Elem(xPath));

}

std::string XmlParameters::GetSubFileName(TiXmlNode* vNode,std::string xPath)
{
	return SubFileName(Elem(vNode,xPath));
}


bool XmlParameters::Exists(std::string xPath)
{
	unsigned number=Numbers(xPath);
	if(number>0)
	{
		return true;
	}

	return false;
}

bool XmlParameters::Exists(TiXmlNode* vNode,std::string xPath)
{
	unsigned number=Numbers(vNode,xPath);
	if(number>0)
	{
		return true;
	}

	return false;
}


void XmlParameters::ExistsOrDie(std::string xPath,std::string vMsg)
{
	if(!Exists(xPath))
	{

//		cout<<"Your parameter "<<xpath<<" is missing: program stopped"<<endl;
//		cout<<msg<<endl;
//		exit(1);
		Error err("exists_or_die",xPath,"Your parameter "+xPath+" is missing. You are in the exist or die function -> die. (Paradox: how could you die if you do not exist? ) \n"+vMsg);
		throw err;
	}
}

void XmlParameters::ExistsOrDie(TiXmlNode* vNode,std::string xPath,std::string vMsg)
{
	if(!Exists(vNode,xPath))
	{
	//	cout<<"Your parameter "<<xpath<<" is missing: program stopped"<<endl;
	//	cout<<msg<<endl;
		//exit(1);

		Error err("exists_or_die",xPath,"Your parameter "+xPath+" is missing. You are in the exist or die function -> die. (Paradox: how could you die if you do not exist?)\n You are in a node-type function, if the error is not in the xml, please check that your xpath starts by // and not by /, it could be the solution.\n"+vMsg);
		throw err;
	}
}


std::string XmlParameters::Process(std::string xPath)
{
	return TinyXPath::S_xpath_string(mpDoc->RootElement(),xPath.c_str());
}


std::string XmlParameters::Process(TiXmlNode* vNode, std::string xPath)
{
	return TinyXPath::S_xpath_string(vNode,xPath.c_str());
}


unsigned XmlParameters::Numbers(std::string xPath)
{
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	return xp_proc.u_compute_xpath_node_set();
}

unsigned XmlParameters::Numbers(TiXmlNode* vNode,std::string xPath)
{
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	return xp_proc.u_compute_xpath_node_set();
}

std::string XmlParameters::Elem(std::string xPath,unsigned vNumber)
{
	xPath+="/text()";
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		Log::mE<<"For the xpath :"<<xPath<<" there are only "<<num<<" elements that are corresponding to your request in the document '"<<mFilename<<"', while you are asking for the "<<vNumber<<" th element"<<endl;
		//exit(1);

		Error err("elem",xPath,"Element not found, check if it exists (and if there are enough  elements");
		throw err;
	}
	return xp_proc.XNp_get_xpath_node(vNumber)->Value();

}

std::string XmlParameters::Elem(TiXmlNode* vNode,std::string xPath,unsigned vNumber)
{
	xPath+="/text()";
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		Log::mE<<"For the xpath :"<<xPath<<" there are only "<<num<<" values that are checking it in you document '"<<mFilename<<"', while you are asking for the "<<vNumber<<" th element"<<endl;
		//exit(1);
		Error err("elem",xPath,"Element not found, check if it exists (and if there are enough  elements. \n You are in a node-type function, if the error is not in the xml, please check that your xpath starts by // and not by /, it could be the solution.");
		throw err;
	}
	return xp_proc.XNp_get_xpath_node(vNumber)->Value();

}


TiXmlNode* XmlParameters::GetNode(std::string xPath,unsigned vNumber)
{
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		Log::mE<<"For the xpath :"<<xPath<<" there are only "<<num<<" values that are checing it in you document '"<<mFilename<<"', while you are asking for the "<<vNumber<<" th element"<<endl;
//		exit(1);

		Error err("get_node",xPath,"Element not found, check if it exists (and if there are enough  elements.");
		throw err;
	}
	return xp_proc.XNp_get_xpath_node(vNumber);
}

TiXmlNode* XmlParameters::GetNode(TiXmlNode* vNode,std::string xPath,unsigned vNumber)
{
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		Log::mE<<"For the xpath :"<<xPath<<" there are only "<<num<<" values that are checing it in you document '"<<mFilename<<"', while you are asking for the "<<vNumber<<" th element"<<endl;
//		exit(1);
		Error err("elem",xPath,"Element not found, check if it exists (and if there are enough  elements. \n You are in a node-type function, if the error is not in the xml, please check that your xpath starts by // and not by /, it could be the solution.");
		throw err;
	}
	return xp_proc.XNp_get_xpath_node(vNumber);
}

std::vector<TiXmlNode*> XmlParameters::GetNodes(TiXmlNode* vNode,std::string xPath)
{
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	vector<TiXmlNode*> resu;
	unsigned num=xp_proc.u_compute_xpath_node_set();
/*	if(num<number)
	{
		cout<<"For the xpath :"<<xpath<<" there are only "<<num<<" values that are checing it in you document '"<<filename<<"', while you are asking for the "<<number<<" th element"<<endl;
		exit(1);
	}
	return xp_proc.XNp_get_xpath_node(number);
	*/
	for(unsigned i=0;i<num;i++)
	{
		resu.push_back(xp_proc.XNp_get_xpath_node(i));
	}
	return resu;
}


std::vector<TiXmlNode*> XmlParameters::GetNodes(std::string xPath)
{
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	vector<TiXmlNode*> resu;
	unsigned num=xp_proc.u_compute_xpath_node_set();
	for(size_t i=0;i<num;i++)
	{
		resu.push_back(xp_proc.XNp_get_xpath_node(i));
	}
	return resu;
}
unsigned XmlParameters::NbParams(std::string xPath,std::string vParam,unsigned vNumberatt)
{
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
//	if(vNumberatt>num)
	if( (num>0&&num<vNumberatt)||num<1)
	{
		return 0;
	}
	return num;
}

unsigned XmlParameters::NbParams(TiXmlNode* vNode,std::string xPath,std::string vParam,unsigned vNumberatt)
{
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
//	if(vNumberatt>num)
	if( (num>0&&num<vNumberatt)||num<1)
	{
		return 0;
	}
	return num;
}

std::string XmlParameters::GetKey(TiXmlNode* vNode,std::string xPath,std::string vParam,unsigned vNumber)
{	
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		//cout<<"For the xpath :"<<xpath<<" there are only "<<num<<" values that are checing it in you document '"<<filename<<"', while you are asking for the "<<number<<" th element"<<endl;
		//exit(1);
		Error err("get_key",xPath,"Element not found");
		throw err;
	}
//	cout<<"valeur : "<<num<<endl;
	return xp_proc.XAp_get_xpath_attribute(vNumber)->Value();
}


std::string XmlParameters::GetKey(std::string xPath,std::string vParam,unsigned vNumber)
{	
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if( (num>0&&num<vNumber)||num<1)
	{
		//cout<<"For the xpath :"<<xpath<<" there are only "<<num<<" values that are checing it in you document '"<<filename<<"', while you are asking for the "<<number<<" th element"<<endl;
		Error err("get_key",xPath,"Element not found");
		throw err;
	}
//	cout<<"valeur : "<<num<<endl;
	return xp_proc.XAp_get_xpath_attribute(vNumber)->Value();

}


bool XmlParameters::KeyExists(std::string xPath,std::string vParam)
{
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(mpDoc->RootElement(),xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if(num<1)
	{
		return false;
	}
	return true;
}

bool XmlParameters::KeyExists(TiXmlNode* vNode,std::string xPath,std::string vParam)
{
	xPath+="/attribute::"+vParam;
	TinyXPath::xpath_processor xp_proc(vNode,xPath.c_str());
	unsigned num=xp_proc.u_compute_xpath_node_set();
	if(num<1)
	{
		return false;
	}
	return true;
}


std::vector<std::string> XmlParameters::GetParamArray(std::string xPath,unsigned vNumber)
{	
	vector<string> result;
	std::string value=Elem(xPath,vNumber);
	if(!MathString::LitToutString(value,result))
	{
		Error err("get_param_array",xPath,"Error while reading table");
		throw err;
	}
	return result;
}
std::deque<std::string> XmlParameters::GetParamDeQue(std::string xPath,unsigned vNumber)
{	
	vector<string> result;
	std::string value=Elem(xPath,vNumber);
	if(!MathString::LitToutString(value,result))
	{
		Error err("get_param_array",xPath,"Error while reading table");
		throw err;
	}
	deque<string> resu;
	for(size_t i=0;i<result.size();++i)
	{
		resu.push_back(result[i]);
	}
	return resu;
}
std::vector<std::string> XmlParameters::GetParamArray(TiXmlNode* vNode,std::string xPath,unsigned vNumber)
{

	vector<string> result;
	std::string value=Elem(vNode,xPath,vNumber);
	if(!MathString::LitToutString(value,result))
	{
		Error err("get_param_array",xPath,"Error while reading table");
		throw err;
	}
	return result;
}

std::deque<std::string> XmlParameters::GetParamDeQue(TiXmlNode* vNode,std::string xPath,unsigned vNumber)
{

	vector<string> result;
	std::string value=Elem(vNode,xPath,vNumber);
	if(!MathString::LitToutString(value,result))
	{
		Error err("get_param_array",xPath,"Error while reading table");
		throw err;
	}
	deque<string> resu;
	//std::copy(result.begin(),result.end(),resu.begin());
	for(size_t i=0;i<result.size();++i)
	{
		resu.push_back(result[i]);
	}
	return resu;
}



void XmlParameters::SetMonteCarloActive()
{
	//Log::SetPriority(Log::CONFIG,"SetMonteCarloActive()");
	Log::mI<<"Monte Carlo processing activated for "<<mFilename<<endl;
	mMonteCarloActivated=true;
}
void XmlParameters::SetMonteCarloInactive()
{
	//Log::SetPriority(Log::CONFIG,"SetMonteCarloInactive()");
	Log::mI<<"Monte Carlo processing disabled for "<<mFilename<<endl;
	mMonteCarloActivated=false;
}

bool XmlParameters::GetMCParams(TiXmlNode* vNode,std::string xPath,double & rPercent,double&rRange,unsigned vNumber,std::string vUncertaintyKeyName)
{
	rPercent=0;
	rRange=0;

	/// Check if globally activated
	if(!mMonteCarloActivated)
	{
		return false;
	}
	/// Check if locally activated
	if(!NbParams(vNode,xPath,"setMC",vNumber)>0)
	{
		return false;
	}
	if(GetKey(vNode,xPath,"setMC",vNumber)!="active")
	{
		return false;
	}
	if(!NbParams(vNode,xPath,vUncertaintyKeyName,vNumber)>0)
	{
		if(vUncertaintyKeyName=="uncertainty")
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"The setMC is active for "<<xPath<<" but the uncertainty is not defined... The MC is not taken into account for the key: "<<vUncertaintyKeyName<<" here."<<endl;
			Log::mW<<"This message can be ignored. And should be ignored if another key is defined: typically if only one of \"uncertainty\" and \"fact_uncertainty\" keys are defined"<<endl;
		}
		return false;
	}

	string el=GetKey(vNode,xPath,vUncertaintyKeyName,vNumber);
	double val=0;
	strton(el,val);
	//Log::SetPriority(Log::DEBUGG,"GetMCParams");
	//	Log::SetPriority(Log::WARNING,"GetMCParams!!");
	if(el.find("%")!=string::npos)
	{
		//		Log::mL<<"You have a percentage MC value :"<<val<<endl;;
		rPercent=val;
	}else
	{
		//		Log::mL<<"You have a range MC value :"<<val<<endl;;
		rRange=val;
	}
	return true;
}


bool XmlParameters::GetMCParams(std::string xPath,double & rPercent,double&rRange,unsigned vNumber,std::string vUncertaintyKeyName)
{
	rPercent=-1;
	rRange=-1;
	/// Check if globally activated
	if(!mMonteCarloActivated)
	{
		return false;
	}
	/// Check if locally activated
	if(!NbParams(xPath,"setMC",vNumber)>0)
	{
		return false;
	}
	if(GetKey(xPath,"setMC",vNumber)!="active")
	{
		return false;
	}
	if(!NbParams(xPath,vUncertaintyKeyName,vNumber)>0)
	{
		if(vUncertaintyKeyName=="uncertainty")
		{
			//Log::SetPriority(Log::WARNING);
			Log::mW<<"The setMC is active for "<<xPath<<" but the uncertainty is not defined... The MC is not taken into account for the key: "<<vUncertaintyKeyName<<" here."<<endl;
			Log::mW<<"This message can be ignored. And should be ignored if another key is defined: typically if only one of \"uncertainty\" and \"fact_uncertainty\" keys are defined"<<endl;
		}
		return false;
	}

	string el=GetKey(xPath,vUncertaintyKeyName,vNumber);
	double val=0;
	strton(el,val);
	val=nabs(val);
//	Log::SetPriority(Log::WARNING,"GetMCParams!!");
//	Log::mL<<el<<" - "<<el.find("%")<<" "<<el.size()<<endl;
	if(el.find("%")!=string::npos)
	{
//		Log::mL<<"You have a percentage MC value :"<<val<<endl;;
		rPercent=val;
	}else
	{
//		Log::mL<<"You have a range MC value :"<<val<<endl;;
		rRange=val;
	}
	return true;
}
/*
void XmlParameters::ApplyMC(double & rValue, double vPercent,double vRange)
{
	double sigma=vRange;
	if(vPercent>0)
	{// We are in the percent case
		sigma=vPercent*rValue/100.;
	}
	MathRandom::GetNormal(rValue,sigma);
	rValue+=sigma;
}

*/


                /*       _\|/_
                         (o o)
                 +----oOO-{_}-OOo----------------------------------------------------+
                 |                                                                   |
                 |Here, we will work on the automatic pointer to XmlParameter values.|
                 |                                                                   |
                 +------------------------------------------------------------------*/



std::vector< std::string> XmlParameters::mFileOpenList;
std::map< std::string, XmlParameters*> XmlParameters::mParamMap;
std::map< std::string, unsigned> XmlParameters::mOpenedMap;


XmlParameters* XmlParameters::AttachFileParameter(std::string vFilename)
{ // This is a static function
	if(find(XmlParameters::mFileOpenList.begin(),XmlParameters::mFileOpenList.end(),vFilename)!=XmlParameters::mFileOpenList.end())
	{// Found
//		cout<<"Pas besoin de lire!!!"<<endl;
		XmlParameters::mOpenedMap[vFilename]+=1;
		return XmlParameters::mParamMap[vFilename];
	}
//	cout<<"Nouveau fichier"<<endl;
	//Throw an exception if there is a problem here...
	XmlParameters* newparam=new XmlParameters(vFilename);
	
	XmlParameters::mFileOpenList.push_back(vFilename);
	XmlParameters::mParamMap[vFilename]=newparam;
	XmlParameters::mOpenedMap[vFilename]=1;
	return newparam;
}

void XmlParameters::DetachFileParameter()
{
//	cout<<"Detachement"<<endl;
	XmlParameters::mOpenedMap[mFilename]-=1;
	if(XmlParameters::mOpenedMap[mFilename]==0)
	{
//		cout<<"On efface vraiment!!!"<<endl;
		vector<string>::iterator it;
		it=find(XmlParameters::mFileOpenList.begin(),XmlParameters::mFileOpenList.end(),mFilename);

		XmlParameters::mFileOpenList.erase(it,it+1);

		// \~fr Pan dans le pied!!!
		// \~en shot in the feet!!!
		delete XmlParameters::mParamMap[mFilename];
	}
}

bool XmlParameters::DoINeedAGarbageCollector()
{
	if(mFileOpenList.size()!=0)
	{
		//Log::SetPriority(Log::INFO);
		Log::mW<<"You need a garbage collector!!!"<<endl;
		return true;
	}
	return false;
}



