/**
 * \ingroup file_utils 
 * \file logging.cpp
 * \brief implements the static logging class
 * Copyright G Gronoff Sept 2009
 * Last modification : $Id: logging.cpp 1111 2010-08-12 19:43:33Z gronoff $
 *
 */
#include "logging.hpp"

using namespace std;


// We initialize the stream as cerr...
std::ostream* Log::msStream=&std::cerr;

Log Log::mL;

Log Log::mD;

Log Log::mE;

Log Log::mW;

Log Log::mI;

Log Log::mS;


bool Log::mbIsOfstream=false;
// We initialize the msMessageLog

std::string Log::msMessageLog="#\t"+static_cast<std::string>(PACKAGE_STRING)+"\n#\tCopyright G Gronoff 2009\n#\t"+"Report bugs : "+static_cast<std::string>(PACKAGE_BUGREPORT)+"\n#\n#\tLaunched :"+Log::ReturnDate()+"\n";


// By default, the priority is set to error, because if something append before  the initialization... It is an error
// Please do not discuss my choices!!!

Log::Log()
{
	mNextMessagePriority=Log::ERROR;
}

#ifdef DEBUG

Log::Priority Log::msCoutPriority=Log::DEBUGG;

#else

Log::Priority Log::msCoutPriority=Log::CONFIG;

#endif

Log::Priority Log::msFilePriority=Log::CONFIG;

Log::Priority Log::msCerrPriority=Log::ERROR;

// We initialize the priority names
const string Log::msPriorityNames[] =
{
	"DEBUG",
	"CONFIG",
	"INFO",
	"SECTION",
	"WARNING",
	"ERROR"
};



void Log::SetPriority(Log::Priority vPriority,std::string vMessage)
{
	if(vPriority!=mNextMessagePriority)
	{

		mNextMessagePriority=vPriority;

		if(mNextMessagePriority>=Log::msCoutPriority)
		{
			cout<<"**** "<<Log::msPriorityNames[mNextMessagePriority]<<endl<<"---- "<<vMessage<<endl;
		}
		if(mNextMessagePriority>=Log::msFilePriority)
		{
			*Log::msStream<<"***************************************"<<endl;
			*Log::msStream<<"**** "<<Log::msPriorityNames[mNextMessagePriority]<<endl<<"---- "<<vMessage<<endl;
			*Log::msStream<<"---------------------------------------"<<endl;
		}


		if(mNextMessagePriority>=Log::msCerrPriority)
		{
			cerr<<"***************************************"<<endl;
			cerr<<"**** "<<Log::msPriorityNames[mNextMessagePriority]<<endl<<"---- "<<vMessage<<endl;
			cerr<<"---------------------------------------"<<endl;
		}
	}
}


void Log::AddVersionInfo(std::string vMessage)
{
	Log::msMessageLog+=vMessage+"\n";
}


void Log::Init(std::string filename,Log::Priority vCoutPr,Log::Priority vFilePr,Log::Priority vCerrPr)
{
	if(filename!="")
	{
		ofstream *of=new ofstream(filename.c_str());
		Log::msStream=of;
	//	*Log::msStream<<"Hello log file"<<endl;
		*Log::msStream<<Log::msMessageLog<<endl;
		*Log::msStream<<endl;
		Log::mbIsOfstream=true;
	}
	Log::msFilePriority=vFilePr;
	Log::msCerrPriority=vCerrPr;
	Log::msCoutPriority=vCoutPr;
	Log::mL.SetPriority(Log::CONFIG);
	Log::mD.SetPriority(Log::DEBUGG);
	Log::mI.SetPriority(Log::INFO);
	Log::mS.SetPriority(Log::SECTION);
	Log::mW.SetPriority(Log::WARNING);
	Log::mE.SetPriority(Log::ERROR);
}




std::string Log::ReturnDate()
{
	time_t timestamp = time(0);
	tm *ltime=localtime(&timestamp);
	std::vector<char> buffer(42);
	strftime(&buffer[0],buffer.size(),"%c",ltime);
	string date(&buffer[0]);
	return date;
}

void Log::Close()
{
	//msStream->close();
	if(Log::mbIsOfstream)
		delete msStream;
	Log::mbIsOfstream=false;
}


