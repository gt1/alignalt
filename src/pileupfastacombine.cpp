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

struct PileupLine
{
	int64_t refpos;
	int64_t refdif;
	std::vector< std::pair<char,uint64_t> > pilefreq;

	PileupLine(std::string const & s)
	{
		std::deque<std::string> const tokens = libmaus2::util::stringFunctions::tokenize(s,std::string("\t"));

		if ( tokens.size() >= 2 )
		{
			refpos = atol(tokens[0].c_str());
			refdif = atol(tokens[1].c_str());

			for ( uint64_t i = 2; i < tokens.size(); ++i )
			{
				std::string const & sp = tokens[i];

				std::deque<std::string> const subtokens = libmaus2::util::stringFunctions::tokenize(sp,std::string(","));
				if ( subtokens.size() == 2 && subtokens[0].size() == 1 )
				{
					std::string const & base = subtokens[0];
					std::string const & snum = subtokens[1];
					int64_t const num = atol(snum.c_str());
					pilefreq.push_back(std::pair<char,uint64_t>(base[0],num));
				}
			}
		}
	}

	uint64_t getDepth() const
	{
		uint64_t depth = 0;
		for ( uint64_t i = 0; i < pilefreq.size(); ++i )
			depth += pilefreq[i].second;
		return depth;
	}
};

std::ostream & operator<<(std::ostream & out, PileupLine const & PL)
{
	out << "PileupLine(" << PL.refpos << "," << PL.refdif;
	for ( uint64_t i = 0; i < PL.pilefreq.size(); ++i )
		out << ",(" << PL.pilefreq[i].first << "," << PL.pilefreq[i].second << ")";
	out << ")";
	return out;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);
		std::string const multifn = arg[0];
		std::string const pileAfn = arg[1];
		std::string const pileBfn = arg[2];
		bool const filter = arg.argPresent("f");
		int64_t const freqfilter =  arg.uniqueArgPresent("F") ? atol(arg["F"].c_str()) : 0;

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

				PileupLine pilea(linea);
				PileupLine pileb(lineb);

				#if 0
				std::cerr << pilea << std::endl;
				std::cerr << pileb << std::endl;
				#endif

				int64_t const maxdepth = std::max(pilea.getDepth(),pileb.getDepth());

				bool const same =
					pilea.pilefreq.size() == 1
					&&
					pileb.pilefreq.size() == 1
					&&
					pilea.pilefreq[0].first == pileb.pilefreq[0].first;

				if ( ((! same) || (! filter)) && (maxdepth >= freqfilter) )
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
