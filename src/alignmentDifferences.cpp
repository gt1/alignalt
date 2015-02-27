#include <libmaus/fastx/StreamFastAReader.hpp>

int main()
{
	try
	{
		libmaus::fastx::StreamFastAReaderWrapper SFQR(std::cin);
		std::vector<std::string> patterns;
		std::vector<std::string> names;
		libmaus::fastx::StreamFastAReaderWrapper::pattern_type pat;

		while ( SFQR.getNextPatternUnlocked(pat) )
		{
			patterns.push_back(pat.spattern);
			names.push_back(pat.sid);
		}
		
		for ( uint64_t i = 1; i < patterns.size(); ++i )	
			if ( patterns[i].size() != patterns[0].size() )
			{
				std::cerr << "The two sequences do not have the same length." << std::endl;
				return EXIT_FAILURE;
			}
			
		if ( ! patterns.size() )
		{
			std::cerr << "no patterns found." << std::endl;
			return EXIT_SUCCESS;
		}
		
		std::vector<uint64_t> P(patterns.size(),0);
		std::set<uint64_t> D;

		for ( uint64_t j = 0; j < patterns[0].size(); ++j )
		{
			bool keep = true;
			
			for ( uint64_t z = 0; z < patterns.size(); ++z )
			{
				bool isbase =
				(
					tolower(patterns[z][j]) == 'a' ||
					tolower(patterns[z][j]) == 'c' ||
					tolower(patterns[z][j]) == 'g' ||
					tolower(patterns[z][j]) == 't' ||
					tolower(patterns[z][j]) == 'n'
				);
				
				if ( isbase )
				{
				
				}
				else
				{
					assert ( patterns[z][j] == '-' );
					keep = false;
				}
			}
						
			if ( keep )
			{
				char const ref = tolower(patterns[0][j]);
				bool dif = false;
				for ( uint64_t z = 1; z < patterns.size(); ++z )
					if ( tolower(patterns[z][j]) != ref )
						dif = true;
				if ( dif )
					D.insert(P[0]);
			
				for ( uint64_t z = 0; z < patterns.size(); ++z )
					patterns[z][P[z]++] = toupper(patterns[z][j]);
			}
		}		

		for ( uint64_t z = 0; z < patterns.size(); ++z )
			patterns[z] = patterns[z].substr(0,P[z]);
			
		std::cout 
			<< names[0] << "\t";
		for ( std::set<uint64_t>::const_iterator ita = D.begin(); ita != D.end(); ++ita )
			std::cout.put(patterns[0][*ita]);
		std::cout.put('\t');
		std::cout
			<< patterns[0] << "\n";
		
		for ( uint64_t z = 1; z < patterns.size(); ++z )
		{
			std::cout << names[z] << "\t";
			std::string const & s0 = patterns[0];
			std::string const & so = patterns[z];
			assert ( s0.size() == so.size() );

			uint64_t ncommon = 0;
			uint64_t nlocal = 0;

			for ( std::set<uint64_t>::const_iterator ita = D.begin(); ita != D.end(); ++ita )
			{
				uint64_t const p = *ita;
			
				// difference in other sequences	
				if ( s0[p] == so[p] )
					std::cout.put('.');
				// difference between ref and this
				else
					std::cout.put(so[p]);
			
				// is this difference common in all regions?
				bool common = true;
				for ( uint64_t j = 2; j < patterns.size(); ++j )
					if ( patterns[j][p] != patterns[1][p] )
						common = false;
						
				if ( common )
					ncommon++;
				else
				{
					if ( s0[p] != so[p] )
						nlocal++;
				}
			}
			std::cout.put('\t');

			uint64_t mat = 0;
			uint64_t mis = 0;

			for ( uint64_t i = 0; i < s0.size(); ++i )
				if ( so[i] == s0[i] )
					++mat;
				else
					++mis;
					
			std::cout << (mat+mis) << "/" << mat << "/" << mis << "/" << static_cast<double>(mat)/(mat+mis) << "/" << ncommon << "/" << nlocal << "\t";
			
			for ( uint64_t i = 0; i < s0.size(); ++i )
				if ( so[i] == s0[i] )
					std::cout.put('.');
				else
					std::cout.put(so[i]);
			
			std::cout.put('\n');
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
