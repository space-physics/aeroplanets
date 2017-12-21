/**
 * \file lines.cpp
 * \brief Implementation of the class to add lines to solar flux grids
 * Copyright G Gronoff Sept 2009
 * Last Modification $Id: lines.cpp 1111 2010-08-12 19:43:33Z gronoff $
 *
 */




#include "lines.hpp"

using namespace std;
using namespace boost::assign;


SolarGridLines::SolarGridLines(XmlParameters* vpParam):mpParameters(vpParam)
{
}


void SolarGridLines::InitStandard()
{
	vector<double> MinLines,MaxLines;
//	MinLines+= 48.3846391,40.8870565,40.8238126,33.6949512,26.6606599,22.3748827,21.2281815,19.6984697,17.6386701,16.2139077,16.103275,15.7169271,12.7000371,12.0975288,12.0250207;
//	MaxLines+=48.3746391, 40.8770565, 40.8138126, 33.6849512, 26.6506599, 22.3648827, 21.2181815, 19.6884697, 17.6286702, 16.2039077, 16.093275, 15.7069271, 12.6900371, 12.0875288, 12.0150207;

	MinLines+= 48.3846390948,43.6433626606,40.8870564769,40.8238126276,33.6949512321,26.6606599028,22.3748826596,21.228181507,20.343278667,19.6984696616,17.6386701455,16.2139077305,16.1032750094,15.7169271308,12.7000370514,12.0975287603,12.0250206898;

	MaxLines+=48.3646390948,43.6233626606,40.8670564769,40.8038126276,33.6749512321,26.6406599028,22.3548826596,21.208181507,20.323278667,19.6784696616,17.6186701455,16.1939077305,16.0832750094,15.6969271308,12.6800370514,12.0775287603,12.0050206898;
	mMinLines.resize(MinLines.size());
	mMaxLines.resize(MinLines.size());
	std::copy(MinLines.begin(),MinLines.end(),mMinLines.begin());
	std::copy(MaxLines.begin(),MaxLines.end(),mMaxLines.begin());
}

void SolarGridLines::InitStandardExtended()
{
	vector<double> MinLines,MaxLines;
	MinLines+=48.3756390948,43.6343626606,40.8780564769,40.8148126276,33.6859512321,26.6516599028,22.3658826596,21.219181507,20.334278667,19.6894696616,17.6296701455,16.2049077305,16.0942750094,15.7079271308,12.6910370514,12.0885287603,12.0160206898,10.1998368554,8.89666352887,8.83952662945,8.00933231063,7.99600893105,7.94361370916,7.48254718803;
	MaxLines+=48.3736390948,43.6323626606,40.8760564769,40.8128126276,33.6839512321,26.6496599028,22.3638826596,21.217181507,20.332278667,19.6874696616,17.6276701455,16.2029077305,16.0922750094,15.7059271308,12.6890370514,12.0865287603,12.0140206898,10.1978368554,8.89466352887,8.83752662945,8.00733231063,7.99400893105,7.94161370916,7.48054718803;

	mMinLines.resize(MinLines.size());
	mMaxLines.resize(MinLines.size());
	std::copy(MinLines.begin(),MinLines.end(),mMinLines.begin());
	std::copy(MaxLines.begin(),MaxLines.end(),mMaxLines.begin());
}
void SolarGridLines::InitUserDefined()
{
	mpParameters->ExistsOrDie("/aero_main/sun/model/lines/Egrid","You should define your lines");
	mpParameters->ExistsOrDie("/aero_main/sun/model/lines/width","You should define your lines width (one width for all the lines)");
	double linewidth=0;
	mpParameters->GetValue("/aero_main/sun/model/lines/width",linewidth);
	ublas::vector<double> meanlines;
	mpParameters->Get1DArray("/aero_main/sun/model/lines/Egrid",meanlines);

	unsigned esize=meanlines.size();
	mMinLines.resize(esize);
	mMaxLines.resize(esize);
	for(unsigned i=0;i<esize;++i)
	{
		mMinLines[i]=meanlines[i]+linewidth/2.;
		mMaxLines[i]=meanlines[i]-linewidth/2.;
	}
	
}


void SolarGridLines::AddSolarLines(ublas::vector<double>& rPhGrMin,ublas::vector<double>& rPhGrMax)
{

	// Only if we have a standard grid for the sun.
	mpParameters->ExistsOrDie("/aero_main/sun/grid/st_grid/use_lines","You must define wether you want or not to use lines in your model");
	unsigned type=0;

	mpParameters->GetNKey("/aero_main/sun/grid/st_grid/use_lines","type",type);

	switch(type)
	{
		case 0:
			Log::mD<<"You do not use extra lines for your solar grid"<<endl;
			return; // Ok, we finish here
			break;
		case 1:
			Log::mD<<"You use the standard lines"<<endl;
			InitStandard();
			break;
		case 2:
			Log::mD<<"You use the  user defined lines"<<endl;
			InitUserDefined();
			break;
		case 3: 
			Log::mD<<"You use the extended standard lines (for the Trans* model"<<endl;
			InitStandardExtended();
			break;
		default:
			Log::mE<<"Lines not defined"<<endl;
			Error err("AddSolarLines","switch(types)","Your line type is not recognized -> pan");
			throw err;
	}

	assert(rPhGrMin.size()==rPhGrMax.size());
	assert(rPhGrMin.size()>0);


	if(rPhGrMin[0]<rPhGrMax[0])
	{
		MathGrid::AddGridLines(mMinLines,mMaxLines,rPhGrMin,rPhGrMax);
	}else
	{
		MathGrid::AddGridLines(mMinLines,mMaxLines,rPhGrMax,rPhGrMin);
	}
}

