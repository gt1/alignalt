/*
    alignalt
    Copyright (C) 2016 German Tischler

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

int main()
{
	try
	{
		libmaus2::fastx::StreamFastAReaderWrapper SFAS(std::cin);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		while ( SFAS.getNextPatternUnlocked(pattern) )
		{
			for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
				pattern.spattern[i] = toupper(pattern.spattern[i]);
			pattern.printMultiLine(std::cout,80);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
