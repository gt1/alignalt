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

std::string readline(std::istream & in)
{
	std::string line;
	std::getline(in,line);
	assert ( line.size() );
	return line;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);
		std::string const multifn = arg[0];
		std::string const pileAfn = arg[1];
		std::string const pileBfn = arg[2];

		libmaus2::aio::InputStreamInstance multiISI(multifn);
		libmaus2::fastx::StreamFastAReaderWrapper SFARW(multiISI);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		std::vector < std::pair<std::string,std::string> > VF;
		while ( SFARW.getNextPatternUnlocked(pattern) )
		{
			for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
				pattern.spattern[i] = toupper(pattern.spattern[i]);
			VF.push_back(std::pair<std::string,std::string>(pattern.sid,pattern.spattern));
		}
		assert ( VF.size() == 2 );

		libmaus2::aio::InputStreamInstance pileAISI(pileAfn);
		libmaus2::aio::InputStreamInstance pileBISI(pileBfn);

		assert ( VF[0].second.size() == VF[1].second.size() );
		for ( uint64_t i = 0; i < VF[0].second.size(); ++i )
		{
			char const ca = VF[0].second[i];
			char const cb = VF[1].second[i];

			bool const capad = (ca == '-');
			bool const cbpad = (cb == '-');

			assert ( ! (capad && cbpad ) );

			if ( capad )
			{
				std::string const lineb = readline(pileBISI);
				std::cout << "::" << lineb << std::endl;
			}
			else if ( cbpad )
			{
				std::string const linea = readline(pileAISI);
				std::cout << linea << "::" << std::endl;
			}
			else
			{
				std::string const linea = readline(pileAISI);
				std::string const lineb = readline(pileBISI);
				std::cout << linea << "::" << lineb << std::endl;
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
