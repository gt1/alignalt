/*
    alignalt
    Copyright (C) 2014-2015 German Tischler
    Copyright (C) 2014-2015 Genome Research Limited

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
#include <libmaus/fastx/StreamFastAReader.hpp>
#include <libmaus/regex/PosixRegex.hpp>

struct TypeAndCoordinates
{
	bool valid;
	std::string type;
	int64_t low;
	int64_t high;
	std::string sequence;
	
	TypeAndCoordinates() : valid(false) {}
	TypeAndCoordinates(
		std::string rtype,
		int64_t rlow,
		int64_t rhigh
	) : valid(true), type(rtype), low(rlow), high(rhigh) {}
	
	bool operator<(TypeAndCoordinates const & O) const
	{
		if ( valid != O.valid )
			return valid < O.valid;
		else if ( type != O.type )
			return type < O.type;
		else if ( low != O.low )
			return low < O.low;
		else
			return high < O.high;
	}
	
	bool operator!=(TypeAndCoordinates const & O) const
	{
		return
			(
				(*this<O)
				||
				(O<*this)
			);
	}
	
	bool operator==(TypeAndCoordinates const & O) const
	{
		return !operator!=(O);
	}
};

std::ostream & operator<<(std::ostream & out, TypeAndCoordinates const & T)
{
	if ( T.valid )
		return out << "TypeAndCoordinates(" << T.type << "," << T.low << "," << T.high << ")";
	else
		return out << "TypeAndCoordinates()";
}

TypeAndCoordinates getTypeAndCoordinates(std::string s)
{
	std::string::size_type j = s.size();
	for ( std::string::size_type i = 0; i < s.size(); ++i )
		if ( isspace(s[i]) )
		{
			j = i;
			break;
		}
		
	s = s.substr(0,j);

	std::vector < std::string::size_type > V;
	for ( std::string::size_type i = 0; i < s.size(); ++i )
		if ( s[i] == '_' )
			V.push_back(i);
	if ( V.size() < 2 )
		return TypeAndCoordinates();
	
	std::string::size_type i0 = V[V.size()-2]+1;
	std::string::size_type i1 = V[V.size()-1]+1;

	for ( std::string::size_type i = i0; i < i1-1; ++i )
		if ( ! isdigit(s[i]) )
			return TypeAndCoordinates();
	for ( std::string::size_type i = i1; i < s.size(); ++i )
		if ( ! isdigit(s[i]) )
			return TypeAndCoordinates();
	
	std::istringstream istr1(std::string(s.begin()+i0,s.begin()+i1-1));
	std::istringstream istr2(std::string(s.begin()+i1,s.end()));
	
	uint64_t u1;
	istr1 >> u1;
	uint64_t u2;
	istr2 >> u2;
	
	if ( ! istr1 || istr1.peek() != std::istream::traits_type::eof() )
		return TypeAndCoordinates();
	if ( ! istr2 || istr2.peek() != std::istream::traits_type::eof() )
		return TypeAndCoordinates();
		
	return TypeAndCoordinates(std::string(s.begin(),s.begin()+i0-1),u1,u2);
}

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
				std::cerr << "The sequences do not have the same length." << std::endl;
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

		std::set<TypeAndCoordinates> TCS;
		std::set<std::string> typeset;
		std::string refseq;
		for ( uint64_t i = 1; i < names.size(); ++i )
		{
			TypeAndCoordinates TC = getTypeAndCoordinates(names[i]);
			TC.sequence = patterns[i];
			TCS.insert(TC);
			typeset.insert(TC.type);
			refseq = TC.sequence;
		}
		
		std::map<TypeAndCoordinates,bool> sameexists;
		bool allfoundsame = true;
		for ( std::set<TypeAndCoordinates>::const_iterator ita = TCS.begin(); ita != TCS.end(); ++ita )
		{
			bool foundsame = false;
			for ( std::set<TypeAndCoordinates>::const_iterator sita = TCS.begin(); sita != TCS.end(); ++sita )
				if ( sita->type != ita->type && ita->sequence == sita->sequence )
					foundsame = true;
			sameexists[*ita] = foundsame;
			allfoundsame = allfoundsame && foundsame;
		}
		std::map<TypeAndCoordinates,bool> differentexists;
		bool allfounddifferent = true;
		for ( std::set<TypeAndCoordinates>::const_iterator ita = TCS.begin(); ita != TCS.end(); ++ita )
		{
			bool founddiff = false;
			for ( std::set<TypeAndCoordinates>::const_iterator sita = TCS.begin(); sita != TCS.end(); ++sita )
				if ( sita->type != ita->type && ita->sequence != sita->sequence )
					founddiff = true;
			differentexists[*ita] = founddiff;
			allfounddifferent = allfounddifferent && founddiff;
		}
		
		std::ostringstream commentstr;
		if ( TCS.size() == 1 )
		{
			std::string const type = TCS.begin()->type;
			commentstr << "nonovlp:" << type << ":[" << TCS.begin()->low << ":" <<  TCS.begin()->high << "]";
		}
		else if ( TCS.size() > 1 )
		{
			// more than one alignment but same region
			if ( typeset.size() == 1 )
			{
				std::string const type = *(typeset.begin());
				std::vector<TypeAndCoordinates> VTCS(TCS.begin(),TCS.end());
				bool sameseq = true;
				std::string const refseq = VTCS[0].sequence;
				for ( size_t i = 1; i < VTCS.size(); ++i )
					sameseq = sameseq && VTCS[i].sequence == refseq;
				commentstr << "multimapsingle:" << type << ":";
				
				for ( uint64_t i = 0; i < VTCS.size(); ++i )
					commentstr << "[" << VTCS[i].low << "," << VTCS[i].high << "]";
				
				commentstr << ":sameseq=" << (sameseq?"yes":"no");
			}
			else if ( typeset.size() > 1 )
			{
				if ( typeset.size() == 2 )
				{
					assert ( TCS.size() >= 2 );
					
					if ( TCS.size() == 2 )
					{
						std::vector<TypeAndCoordinates> VTCS(TCS.begin(),TCS.end());
						std::string const typea = VTCS[0].type;
						uint64_t const lowa = VTCS[0].low;
						uint64_t const higha = VTCS[0].high;
						std::string const typeb = VTCS[1].type;
						uint64_t const lowb = VTCS[1].low;
						uint64_t const highb = VTCS[1].high;
						bool const sameseq = VTCS[0].sequence == VTCS[1].sequence;
						commentstr << "ovlp:" << typea << ":[" << lowa << "," << higha << "]:" << typeb << ":[" << lowb << "," << highb << "]:" << (sameseq ? "same":"dif");
					}
					else
					{
						std::vector<std::string> Vtypeset(typeset.begin(),typeset.end());
						assert ( TCS.size() > 2 );
						commentstr << "multiovlp:" << Vtypeset[0] << ":";

						std::vector<TypeAndCoordinates> VTCS(TCS.begin(),TCS.end());
						for ( uint64_t i = 0; i < VTCS.size(); ++i )
							if ( VTCS[i].type == Vtypeset[0] )
								commentstr << "[" << VTCS[i].low << "," << VTCS[i].high << "]";
						
						commentstr << ":" << Vtypeset[1] << ":";
						for ( uint64_t i = 0; i < VTCS.size(); ++i )
							if ( VTCS[i].type == Vtypeset[1] )
								commentstr << "[" << VTCS[i].low << "," << VTCS[i].high << "]";
						commentstr << ":allfoundsame=" << (allfoundsame?"yes":"no") << ":allfounddiff=" << (allfounddifferent?"yes":"no");
					}
				}
			}
		}
	
		bool sameseq = names.size()>1;
		for ( std::set<TypeAndCoordinates>::const_iterator ita = TCS.begin(); ita != TCS.end(); ++ita )
		{
			sameseq = sameseq && ita->sequence == refseq;
			// std::cerr << *ita << "\tsameexists=" << sameexists[*ita] << "\tdifferentexists=" << differentexists[*ita] << std::endl;
		}
					
		std::cout 
			<< names[0] << "\t";
		for ( std::set<uint64_t>::const_iterator ita = D.begin(); ita != D.end(); ++ita )
			std::cout.put(patterns[0][*ita]);
		std::cout.put('\t');
		std::cout
			<< patterns[0] << "\t" << commentstr.str() << "\n";
		
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
