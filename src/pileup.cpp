/*
    alignalt
    Copyright (C) 2009-2016 German Tischler
    Copyright (C) 2011-2015 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/util/ArgInfo.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type Pdec(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
		libmaus2::bambam::BamAlignmentDecoder & decoder = Pdec->getDecoder();
		libmaus2::bambam::BamHeader const & header = decoder.getHeader();		
		std::vector< libmaus2::bambam::PileVectorElement > GPV;

		assert ( header.getNumRef() >= 1 );

		std::string const refname = header.getRefIDName(0);
		uint64_t const reflen = header.getRefIDLength(0);

		while ( decoder.readAlignment() )
		{
			libmaus2::bambam::BamAlignment const & algn = decoder.getAlignment();
			if ( algn.isMapped() )
			{
				std::vector< libmaus2::bambam::PileVectorElement > const PV = algn.getPileVector();
				for ( uint64_t i = 0; i < PV.size(); ++i )
					GPV.push_back(PV[i]);
			}
		}
		
		std::sort(GPV.begin(),GPV.end());
		
		uint64_t ilow = 0;
		int64_t const minborderdist = 10;
		uint64_t const mindepth = 5;
		uint64_t next = 0;
		
		std::cout << ">" << refname << "\n";
		
		while ( ilow < GPV.size() )
		{
			uint64_t ihigh = ilow+1;
			
			while ( 
				ihigh < GPV.size() &&
				GPV[ilow].refid == GPV[ihigh].refid &&
				GPV[ilow].refpos == GPV[ihigh].refpos &&
				GPV[ilow].predif == GPV[ihigh].predif 
			)
				++ihigh;
			
			std::vector < libmaus2::bambam::PileVectorElement > GPVfilt;
			for ( uint64_t i = ilow; i < ihigh; ++i )
				if ( GPV[i].readpos >= minborderdist && GPV[i].readbackpos >= minborderdist )
					GPVfilt.push_back(GPV[i]);
			
			if ( GPVfilt.size() >= mindepth && GPVfilt.size() )
			{
				std::map<char,uint64_t> S;
				
				for ( uint64_t i = 0; i < GPVfilt.size(); ++i )
					S [ GPVfilt[i].sym ] ++;
					
				std::vector < std::pair<uint64_t,char> > SV;
				for ( std::map<char,uint64_t>::const_iterator ita = S.begin(); ita != S.end(); ++ita )
					SV.push_back(std::pair<uint64_t,char>(ita->second,ita->first) );
				
				std::sort(SV.begin(),SV.end());
				std::reverse(SV.begin(),SV.end());
				
				#if 0
				std::cerr << std::string(80,'-') << std::endl;
				for ( uint64_t i = 0; i < GPVfilt.size(); ++i )
					std::cerr << GPVfilt[i] << std::endl;
				for ( uint64_t i = 0; i < SV.size(); ++i )
					std::cerr << SV[i].second << " " << SV[i].first << std::endl;
				#endif
				
				if ( SV.size() == 1 && GPVfilt[0].predif == 0 )
				{
					while ( next < GPVfilt[0].refpos )
					{
						std::cout.put('N');
						next += 1;
					}
					
					std::cerr << GPVfilt[0].refpos << "\t" << GPVfilt[0].predif << "\t" << SV[0].second << "\t" << SV[0].first << std::endl;
					
					std::cout.put(SV[0].second);
					next += 1;
				}
			}

			ilow = ihigh;
		}

		while ( next < reflen )
		{
			std::cout.put('N');
			next += 1;
		}

		std::cout.put('\n');
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
