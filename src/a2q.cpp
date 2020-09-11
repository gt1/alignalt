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
#include <iostream>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/fastx/FastQReader.hpp>

int main(int argc, char *argv[])
{
	try
	{
		::libmaus2::util::ArgInfo arginfo(argc,argv);
		typedef ::libmaus2::fastx::FastAReader reader_type;
		typedef ::libmaus2::fastx::FastQReader q_reader_type;
		typedef reader_type::unique_ptr_type reader_ptr_type;
		typedef reader_type::pattern_type pattern_type;
		typedef q_reader_type::pattern_type q_pattern_type;
		reader_ptr_type reader;

		if ( arginfo.restargs.size() < 1 )
		{
			::libmaus2::exception::LibMausException se;
			se.getStream() << "usage: " << arginfo.progname << " <inputfile.fa>";
			se.finish();
			throw se;
		}

		if ( arginfo.restargs[0] == "-" )
			reader = std::move(reader_ptr_type(new reader_type("/dev/stdin")));
		else
			reader = std::move(reader_ptr_type(new reader_type(arginfo.restargs[0])));

		pattern_type pattern;
		q_pattern_type q_pattern;

		std::cerr << "Parsing file...";
		while ( reader->getNextPatternUnlocked(pattern) )
		{
			q_pattern.spattern = ::libmaus2::fastx::remapString(::libmaus2::fastx::mapString(pattern.spattern));
			q_pattern.pattern = q_pattern.spattern.c_str();
			q_pattern.sid = pattern.sid;

			for ( uint64_t i = 0; i < q_pattern.sid.size(); ++i )
				if ( isspace(q_pattern.sid[i]) )
					q_pattern.sid[i] = '_';
				#if 0
				else if ( q_pattern.sid[i] == '=' )
					q_pattern.sid[i] = '_';
				#endif

			q_pattern.plus = "";
			q_pattern.quality = std::string(pattern.spattern.size(),'I');

			std::cout << q_pattern;

			if ( (pattern.getPatID() & (1024*1024-1)) == 0 )
				std::cerr << "(" << (pattern.getPatID()/(1024*1024)) << "m)";
		}
		std::cerr << "done." << std::endl;

	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}
}
