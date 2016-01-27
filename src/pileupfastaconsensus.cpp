/*
    alignalt
    Copyright (C) 2009-2016 German Tischler

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
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/stringFunctions.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);
		std::string const infofn = arg[0];

		libmaus2::aio::InputStreamInstance ISI(infofn);
		std::vector< std::vector<uint64_t> > VINTV;
		while ( ISI )
		{
			std::string line;
			std::getline(ISI,line);
			if ( line.size() )
			{
				std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(line,std::string("\t"));
				assert ( tokens.size() >= 8 );
				if ( tokens.size() >= 8 )
				{
					std::string clen = tokens[7];
					std::deque<std::string> ctokens = libmaus2::util::stringFunctions::tokenize(clen,std::string(";"));
					std::vector<uint64_t> VV;
					for ( uint64_t i = 0; i < ctokens.size(); ++i )
						VV.push_back(atol(ctokens[i].c_str()));
					VINTV.push_back(VV);
				}
			}
		}


		libmaus2::fastx::StreamFastAReaderWrapper SFARW(std::cin);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		std::vector < std::pair<std::string,std::string> > VF;
		uint64_t cc = 0;
		while ( SFARW.getNextPatternUnlocked(pattern) )
		{
			for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
				pattern.spattern[i] = toupper(pattern.spattern[i]);
			if ( cc++ )
				VF.push_back(std::pair<std::string,std::string>(pattern.sid,pattern.spattern));
		}
		assert ( VF.size() );
		assert ( VF.size() == VINTV.size() );
		for ( uint64_t i = 1; i < VF.size(); ++i )
			assert ( VF[i].second.size() == VF[0].second.size() );

		std::set < std::pair<uint64_t,uint64_t> > active;
		for ( uint64_t i = 0; i < VINTV.size(); ++i )
		{
			std::string const & line = VF[i].second;
			uint64_t p = 0;

			for ( uint64_t j = 0; j < VINTV[i].size(); ++j )
			{
				uint64_t l = VINTV[i][j];

				while ( p < line.size() && line[p] == '-' )
					++p;

				while ( p < line.size() && l )
				{
					if ( line[p] != '-' )
					{
						active.insert(std::pair<uint64_t,uint64_t>(i,p));
						--l;
					}
					p += 1;
				}

				assert ( ! l );
			}
		}

		uint64_t maxdepth = 0;
		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
		{
			uint64_t c = 0;
			for ( uint64_t j = 0; j < VF.size(); ++j )
				if ( VF[j].second[i] != '-' )
					++c;
			maxdepth = std::max(maxdepth,c);
		}

		std::ostringstream consensusdatastr;
		std::ostringstream consensusnamestr;
		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
		{
			std::map<char,uint64_t> M;
			for ( uint64_t j = 0; j < VF.size(); ++j )
				if (
					active.find(std::pair<uint64_t,uint64_t>(j,i))
					!= active.end()
				)
				{
					M [ VF[j].second[i] ] ++;
				}

			std::vector < std::pair<uint64_t,char> > VC;
			for ( std::map<char,uint64_t>::const_iterator ita = M.begin(); ita != M.end(); ++ita )
				VC.push_back(std::pair<uint64_t,char>(ita->second,ita->first));
			std::sort(VC.begin(),VC.end());
			std::reverse(VC.begin(),VC.end());

			if ( VC.size() && VC[0].second != '-' && VC[0].first >= maxdepth/2 )
			{
				consensusdatastr.put(VC[0].second);
				//consensusnamestr << "_" << i;
			}
			else
			{
				consensusdatastr.put('N');
			}
		}

		std::cout << ">consensus" << consensusnamestr.str() << std::endl;
		std::cout << consensusdatastr.str() << std::endl;

		#if 0
		uint64_t reflen = 0;
		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
			if ( VF[0].second[i] != '-' )
				reflen += 1;

		int64_t refpos = reflen;
		int64_t refdif = 0;
		std::map<uint64_t, std::pair<uint64_t,int64_t> > PM;
		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
		{
			uint64_t const p = VF[0].second.size()-i-1;

			if ( VF[0].second[p] == '-' )
				refdif -= 1;
			else
			{
				refdif = 0;
				refpos -= 1;
			}

			// std::cerr << p << " " << refpos << " " << refdif << std::endl;
			PM[p] = std::pair<uint64_t,int64_t>(refpos,refdif);
		}

		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
		{
			std::map<char,uint64_t> M;
			for ( uint64_t j = 1; j < VF.size(); ++j )
				if (
					active.find(std::pair<uint64_t,uint64_t>(j,i))
					!= active.end()
				)
					M [ VF[j].second[i] ] ++;
			std::vector < std::pair<uint64_t,char> > VC;
			for ( std::map<char,uint64_t>::const_iterator ita = M.begin(); ita != M.end(); ++ita )
				VC.push_back(std::pair<uint64_t,char>(ita->second,ita->first));
			std::sort(VC.begin(),VC.end());
			std::reverse(VC.begin(),VC.end());

			std::pair<uint64_t,int64_t> PP = PM.find(i)->second;

			std::cout << PP.first << "\t" << PP.second;
			for ( uint64_t i = 0; i < VC.size(); ++i )
				std::cout << "\t" << VC[i].second << "," << VC[i].first;
			std::cout << std::endl;
		}
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
