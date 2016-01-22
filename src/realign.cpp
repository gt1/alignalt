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
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/fastx/FastAStreamSet.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/math/IntegerInterval.hpp>

struct Contig
{
	std::string name;
	std::string data;
	std::string frontclipped;
	std::string backclipped;
	uint64_t low;
	uint64_t high;
	uint64_t clow;
	uint64_t chigh;

	Contig(std::string const & rname = std::string(), std::string const & rdata = std::string(),
		std::string const & rfrontclipped = std::string(),
		std::string const & rbackclipped = std::string(),
		uint64_t const rlow = 0, uint64_t const rhigh = 0,
		uint64_t const rclow = 0, uint64_t const rchigh = 0)
	: name(rname), data(rdata), frontclipped(rfrontclipped), backclipped(rbackclipped), low(rlow), high(rhigh), clow(rclow), chigh(rchigh)
	{

	}
};

std::ostream & operator<<(std::ostream & out, Contig const & contig)
{
	out << "Contig(" << contig.name << "," << contig.data << "," << contig.frontclipped << "," << contig.backclipped << "," << contig.low << "," << contig.high << "," << contig.clow << "," << contig.chigh << ")";
	return out;
}

struct Compactify
{
	template<typename iterator, typename count_type = uint64_t>
	static std::vector< std::pair< typename std::iterator_traits<iterator>::value_type, count_type > > compactify(iterator ita, iterator ite)
	{
		std::vector< std::pair< typename std::iterator_traits<iterator>::value_type, count_type > > V;
		while ( ita != ite )
		{
			uint64_t c = 1;
			typename std::iterator_traits<iterator>::value_type const v = *(ita);
			++ita;

			while ( ita != ite && *ita == v )
			{
				++ita;
				++c;
			}

			V.push_back(std::pair< typename std::iterator_traits<iterator>::value_type, count_type >(v,c));
		}

		return V;
	}
};

std::vector<std::string> tokenize(std::string const & sintv, char const br = ';')
{
	std::vector<std::string> subtokens;
	uint64_t ilow = 0;
	while ( ilow < sintv.size() )
	{
		uint64_t const ihigh = sintv.find(br,ilow);

		if ( ihigh == std::string::npos )
		{
			subtokens.push_back(sintv.substr(ilow));
			ilow = sintv.size();
		}
		else
		{
			subtokens.push_back(sintv.substr(ilow,ihigh-ilow));
			ilow = ihigh+1;
		}
	}

	return subtokens;
}

struct Info
{
	std::string seqid;
	std::string refid;
	double identity;
	uint64_t equal;
	uint64_t diff;
	std::vector<Contig > Vintv;

	void init(std::string const & s)
	{
		std::vector<std::string> tokens;

		uint64_t low = 0;
		while ( low != s.size() )
		{
			while ( low != s.size() && isspace(s[low]) )
				++low;
			if ( low != s.size() )
			{
				uint64_t high = low+1;
				while ( high != s.size() && !isspace(s[high]) )
					++high;

				tokens.push_back(s.substr(low,high-low));

				low = high;
			}
		}

		if ( tokens.size() == 6 )
		{
			#if 0
			for ( uint64_t i = 0; i < tokens.size(); ++i )
				std::cerr << "token " << tokens[i] << std::endl;
			std::cerr << std::endl;
			#endif

			std::string const sintv = tokens[5];

			std::vector<std::string> subtokens = tokenize(sintv,';');

			for ( uint64_t j = 0; j < subtokens.size(); ++j )
			{
				std::string const subt = subtokens[j];
				std::vector<std::string> subsubtokens = tokenize(subt,':');

				#if 0
				for ( uint64_t k = 0; k < subsubtokens.size(); ++k )
					std::cerr << "[" << k << "]=" << subsubtokens[k] << std::endl;
				#endif

				if ( subsubtokens.size() >= 6 )
				{
					std::string const contig = subsubtokens[0];
					std::string const contigdata = subsubtokens[1];
					std::string const frontclipped = subsubtokens[2];
					std::string const backclipped = subsubtokens[3];

					std::istringstream istr0(subsubtokens[4]);
					libmaus2::math::IntegerInterval<int64_t> intv0(istr0);
					//std::cerr << intv0 << std::endl;

					std::istringstream istr1(subsubtokens[5]);
					libmaus2::math::IntegerInterval<int64_t> intv1(istr1);
					//std::cerr << intv1 << std::endl;

					Vintv.push_back(Contig(contig,contigdata,frontclipped,backclipped,intv0.from,intv0.to+1,intv1.from,intv1.to+1));
				}
			}

			diff = atol(tokens[4].c_str());
			equal = atol(tokens[3].c_str());
			identity = atof(tokens[2].c_str());
			refid = tokens[1];
			seqid = tokens[0];

			#if 0
			std::cerr << "seqid=" << seqid << std::endl;
			std::cerr << "refid=" << refid << std::endl;
			std::cerr << "identity=" << identity << std::endl;
			std::cerr << "equal=" << equal << std::endl;
			std::cerr << "diff=" << diff << std::endl;
			for ( uint64_t i = 0; i < Vintv.size(); ++i )
				std::cerr << Vintv[i] << std::endl;
			#endif
		}
	}

	Info()
	{}

	Info(std::string const & s)
	{
		init(s);
	}

	Info(std::istream & in)
	{
		std::string s;
		std::getline(in,s);
		init(s);
	}
};
// ABC8_4_1_000041156800_I24       region_a.fa     0.886454        445     57      [2917,3193];[8165,8389];

#include <libmaus2/bambam/BamFlagBase.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::util::ArgParser const arg(argc,argv);

		libmaus2::aio::InputStreamInstance INFOISI(arg[1]);
		std::vector<Info> Vinfo;
		while ( INFOISI )
		{
			Info info(INFOISI);
			if ( info.seqid.size() )
				Vinfo.push_back(info);
		}

		libmaus2::aio::InputStreamInstance FASS0ISI(arg[0]);
		std::vector< std::pair<std::string,std::string> > Vseq;
		#if 0
		libmaus2::fastx::FastAStreamSet FASS0(FASS0ISI);
		std::pair<std::string,libmaus2::fastx::FastAStream::shared_ptr_type> P;
		while ( FASS0.getNextStream(P) )
		{
			std::istream & in = *(P.second);
			std::string const & name = P.first;
			libmaus2::autoarray::AutoArray<char> B(8192,false);
			std::ostringstream ostr;
			while ( in )
			{
				in.read(B.begin(),B.size());
				ostr.write(B.begin(),in.gcount());
			}

			Vseq.push_back(std::pair<std::string,std::string>(name,ostr.str()));
		}
		#endif
		libmaus2::fastx::StreamFastAReaderWrapper SFASW(FASS0ISI);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		while ( SFASW.getNextPatternUnlocked(pattern) )
		{
			Vseq.push_back(std::pair<std::string,std::string>(pattern.sid,pattern.spattern));
		}

		assert ( Vseq.size() == Vinfo.size() + 1 );
		libmaus2::bambam::BamSeqEncodeTable seqenc;

		std::ostringstream bamheaderstr;

		uint64_t refseqlen = 0;
		std::ostringstream unpaddedrefstr;
		for ( uint64_t i = 0; i < Vseq[0].second.size(); ++i )
			if ( Vseq[0].second[i] != '-' )
			{
				refseqlen += 1;
				unpaddedrefstr.put(Vseq[0].second[i]);
			}
		std::string const unpaddedref = unpaddedrefstr.str();

		bamheaderstr << "@HD\tVN:1.5\tSO:unknown\n";
		bamheaderstr << "@SQ\tSN:" << Vseq[0].first << "\t" << "LN:" << refseqlen << "\n";


		std::string const bamheadertext = bamheaderstr.str();
		libmaus2::bambam::BamHeader const bamheader(bamheadertext);

		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type blockwriter(libmaus2::bambam::BamBlockWriterBaseFactory::construct(bamheader,arginfo));

		// iterate over clones
		for ( uint64_t z = 1; z < Vseq.size(); ++z )
		{
			// reference line
			std::string const ref = Vseq[0].second;
			// clone line
			std::string const oth = Vseq[z].second;
			// info
			Info & info = Vinfo[z-1];
			// check name matches
			assert ( Vseq[z].first == Vinfo[z-1].seqid );

			std::map<uint64_t,uint64_t> contigstart, contigend;

			struct ContigInterval
			{
				uint64_t ciglow;
				uint64_t cighigh;

				uint64_t startrefpos;
				uint64_t startcpos;
			};

			std::map<uint64_t,ContigInterval> contigintv;

			for ( uint64_t i = 0; i < info.Vintv.size(); ++i )
			{
				Contig const & contig = info.Vintv[i];
				contigstart[contig.clow] = i;
				contigend[contig.chigh] = i;
			}

			#if 0
			std::cerr << "---" << std::endl;
			std::cerr << ref.substr(0,80) << std::endl;
			std::cerr << oth.substr(0,80) << std::endl;
			#endif

			assert ( ref.size() == oth.size() );

			std::vector < libmaus2::bambam::BamFlagBase::bam_cigar_ops > cigops;
			uint64_t refpos = 0;
			uint64_t cpos = 0;
			int64_t activecontig = -1;
			std::vector<std::string> recs(info.Vintv.size(),std::string());

			for ( uint64_t i = 0; i < ref.size(); ++i )
			{
				bool const refpad = (ref[i] == '-');
				bool const othpad = (oth[i] == '-');

				if ( ! refpad )
					assert ( ref[i] == unpaddedref[refpos] );

				if ( contigend.find(cpos) != contigend.end() )
				{
					uint64_t const contigid = contigend.find(cpos)->second;
					//std::cerr << "setting end for contig " << contigid << std::endl;
					if ( contigintv[contigid].ciglow == contigintv[contigid].cighigh )
						contigintv[contigid].cighigh = cigops.size();
					activecontig = -1;
				}

				if ( contigstart.find(cpos) != contigstart.end() )
				{
					uint64_t const contigid = contigstart.find(cpos)->second;
					if ( contigintv.find(contigid) == contigintv.end() )
					{
						contigintv[contigid] = ContigInterval();
						contigintv[contigid].ciglow = cigops.size();
						contigintv[contigid].cighigh = cigops.size();
						contigintv[contigid].startrefpos = refpos;
						contigintv[contigid].startcpos = cpos;
					}
					activecontig = contigid;
				}

				if ( (! othpad) && (activecontig != -1) )
					recs[activecontig] += oth[i];

				if ( refpad && othpad )
				{
					// intentionally left blank
				}
				// insert character into reference
				else if ( (!othpad) && refpad )
				{
					cigops.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS);
					cpos += 1;
				}
				// delete character from reference
				else if ( (!refpad) && othpad )
				{
					cigops.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL);
					refpos += 1;
				}
				else if ( ref[i] == oth[i] )
				{
					cigops.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL);
					refpos += 1;
					cpos += 1;
				}
				else
				{
					assert ( ! refpad );
					assert ( ! othpad );
					assert ( ref[i] != oth[i] );
					cigops.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF);
					refpos += 1;
					cpos += 1;
				}
			}

			if ( contigend.find(cpos) != contigend.end() )
			{
				uint64_t const contigid = contigend.find(cpos)->second;
				//std::cerr << "setting end for contig " << contigid << std::endl;
				if ( contigintv[contigid].ciglow == contigintv[contigid].cighigh )
					contigintv[contigid].cighigh = cigops.size();
				activecontig = -1;
			}

			assert ( contigstart.size() == info.Vintv.size() );

			for ( uint64_t c = 0; c < info.Vintv.size(); ++c )
			{
				//std::cerr << "c=" << c << " .size()=" << info.Vintv.size() << std::endl;

				assert ( contigintv.find(c) != contigintv.end() );
				ContigInterval const & CI = contigintv.find(c)->second;
				assert ( recs[c] == info.Vintv[c].data );
				Contig & contig = info.Vintv[c];

				std::vector < std::pair<libmaus2::bambam::BamFlagBase::bam_cigar_ops,uint64_t> > const cigcompact =
					Compactify::compactify<
						std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops>::const_iterator
					>(cigops.begin() + CI.ciglow,cigops.begin() + CI.cighigh);

				uint64_t cc = 0;
				for ( uint64_t i = 0; i < cigcompact.size(); ++i )
					switch ( cigcompact[i].first )
					{
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
							cc += cigcompact[i].second;
							break;
						default:
							break;
					}

				//std::cerr << "cc=" << cc << " expt=" << info.Vintv[c].data.size() << std::endl;
				assert ( cc == info.Vintv[c].data.size() );

				std::string const bamname =
					info.seqid + "_" + contig.name;
				// std::cerr << bamname << std::endl;

				::libmaus2::fastx::EntityBuffer<uint8_t,libmaus2::bambam::BamAlignment::D_array_alloc_type> bambuffer;

				std::string contigdata =
					contig.frontclipped +
					contig.data +
					contig.backclipped;

				std::vector < std::pair<libmaus2::bambam::BamFlagBase::bam_cigar_ops,uint64_t> > cigfinal;
				if ( contig.frontclipped.size() )
					cigfinal.push_back(
						std::pair<libmaus2::bambam::BamFlagBase::bam_cigar_ops,uint64_t>(
							libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,
							contig.frontclipped.size()
						)
					);

				for ( uint64_t i = 0; i < cigcompact.size(); ++i )
					cigfinal.push_back(cigcompact[i]);

				if ( contig.backclipped.size() )
					cigfinal.push_back(
						std::pair<libmaus2::bambam::BamFlagBase::bam_cigar_ops,uint64_t>(
							libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,
							contig.backclipped.size()
						)
					);

				std::string qual(contigdata.size(),255);

				libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment(
					bambuffer,
					seqenc,
					bamname.begin(),
					bamname.size(),
					0 /* ref id */,
					CI.startrefpos,
					255, /* mapping quality */
					0 /* flags */,
					cigfinal.begin(),
					cigfinal.size(),
					-1 /* next ref id */,
					-1 /* next pos */,
					contigdata.size() /* tlen */,
					contigdata.begin(),
					contigdata.size(),
					qual.begin(),
					0, /* qual offset */
					true /* reset buffer */
				);

				libmaus2::bambam::BamAlignment outalgn;
				outalgn.blocksize = bambuffer.swapBuffer(outalgn.D);
				outalgn.checkAlignment();

				blockwriter->writeAlignment(outalgn);
			}
		}

		blockwriter.reset();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
