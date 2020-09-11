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
#include <config.h>

#include <libmaus2/alignment/SimpleLocalAligner.hpp>

#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/bambam/BamFlagBase.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/bambam/CigarOperation.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/StreamFastQReader.hpp>
#include <libmaus2/lcs/MetaLocalEditDistance.hpp>
#include <libmaus2/lcs/MetaEditDistance.hpp>
#include <libmaus2/lf/ImpCompactHuffmanWaveletLF.hpp>
#include <libmaus2/lz/BufferedGzipStream.hpp>
#include <libmaus2/suffixtree/CompressedSuffixTree.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/util/ToUpperTable.hpp>

struct ZSortEntry;
std::ostream & operator<<(std::ostream & out, ZSortEntry const & Z);

struct ZSortEntry
{
	uint64_t zstart;
	uint64_t zend;
	bool zreverse;
	uint64_t seq;
	uint64_t seqstart;
	uint64_t seqend;
	std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops> cigops;

	ZSortEntry()
	{

	}

	ZSortEntry(
		uint64_t const rzstart,
		uint64_t const rzend,
		bool const rzreverse,
		uint64_t const rseq,
		uint64_t const rseqstart,
		uint64_t const rseqend,
		std::vector< libmaus2::bambam::BamFlagBase::bam_cigar_ops  > const & rcigops
	) : zstart(rzstart), zend(rzend), zreverse(rzreverse), seq(rseq), seqstart(rseqstart), seqend(rseqend), cigops(rcigops)
	{
		uint64_t zlen = 0, seqlen = 0;
		for ( uint64_t i = 0; i < cigops.size(); ++i )
			switch ( cigops[i] )
			{
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
					zlen++, seqlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
					zlen++, seqlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
					zlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
					seqlen++;
					break;
				default:
					break;
			}

		if ( zend-zstart != zlen )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "zend-zstart=" << zend-zstart << " != " << "zlen=" << zlen << '\n';
			lme.finish();
			throw lme;
		}
		if (  seqend-seqstart != seqlen )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "seqend-seqstart =" << seqend-seqstart  << " != " << "seqlen=" << seqlen << '\n';
			lme.finish();
			throw lme;
		}
	}

	uint64_t getZStart() const
	{
		return zstart;
	}

	uint64_t getZEnd() const
	{
		return zend;
	}

	std::pair < ZSortEntry, ZSortEntry > zsplit(uint64_t const zmid) const
	{
		// std::cerr << "Splitting " << *this << " at " << zmid << std::endl;

		uint64_t cit = 0;
		uint64_t zlen = 0;
		uint64_t seqlen = 0;

		while ( zlen != zmid-zstart )
		{
			switch ( cigops[cit++] )
			{
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
					zlen++, seqlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
					zlen++, seqlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
					zlen++;
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
					seqlen++;
					break;
				default:
					break;
			}
		}

		if ( !zreverse )
		{
			ZSortEntry left(
				zstart,
				zmid,
				zreverse,
				seq,
				seqstart,
				seqstart+seqlen,
				std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(cigops.begin(),cigops.begin()+cit)
			);
			ZSortEntry right(
				zmid,
				zend,
				zreverse,
				seq,
				seqstart+seqlen,
				seqend,
				std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(cigops.begin()+cit,cigops.end())
			);

			return std::pair < ZSortEntry, ZSortEntry >(left,right);
		}
		else
		{
			ZSortEntry left(
				zstart,
				zmid,
				zreverse,
				seq,
				seqend-seqlen,
				seqend,
				std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(cigops.begin(),cigops.begin()+cit)
			);
			ZSortEntry right(
				zmid,
				zend,
				zreverse,
				seq,
				seqstart,
				seqend-seqlen,
				std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(cigops.begin()+cit,cigops.end())
			);

			return std::pair < ZSortEntry, ZSortEntry >(left,right);
		}
	}

	bool operator<(ZSortEntry const & O) const
	{
		if ( zreverse != O.zreverse )
			return zreverse < O.zreverse;
		else if ( zstart != O.zstart )
			return zstart < O.zstart;
		else if ( zend != O.zend )
			return zend < O.zend;
		else if ( seq != O.seq )
			return seq < O.seq;
		else if ( seqstart != O.seqstart )
			return seqstart < O.seqstart;
		else
			return seqend < O.seqend;
	}
};

std::ostream & operator<<(std::ostream & out, ZSortEntry const & ZSE)
{
	return out << "ZSortEntry("
		<< "zreverse=" << ZSE.zreverse << ","
		<< "zstart=" << ZSE.zstart << ","
		<< "zend=" << ZSE.zend << ","
		<< "seq=" << ZSE.seq << ","
		<< "seqstart=" << ZSE.seqstart << ","
		<< "seqend=" << ZSE.seqend << ")";
}

static uint64_t absdif(uint64_t a, uint64_t b)
{
	if ( a >= b )
		return a-b;
	else
		return b-a;
}

/*
 * wavelet tree based method based on backward search
 */
template<typename _lf_type>
int alignalt(libmaus2::util::ArgInfo const & arginfo, std::string const suffix)
{
	typedef _lf_type lf_type;
	std::string const prefix = arginfo.getUnparsedRestArg(0);
	uint64_t const kmerlen = arginfo.getValue<uint64_t>("kmerlen",16);
	double const errrate = arginfo.getValue<double>("errrate",0.05);
	double const hitfrac = arginfo.getValue<double>("hitfrac",0.8);
	uint64_t const maxkmerfreqthres = arginfo.getValueUnsignedNumeric("maxkmerfreqthres",2048);
	uint64_t const fraglen = arginfo.getValueUnsignedNumeric("fraglen",300);
	uint64_t const clipthres = arginfo.getValueUnsignedNumeric("clipthres",3);
	uint64_t const difthres = arginfo.getValueUnsignedNumeric("difthres",5);

	std::cerr << "[V] constructing aligner...";
	libmaus2::alignment::SimpleLocalAligner<lf_type> aligner(prefix,suffix,kmerlen,errrate,hitfrac,maxkmerfreqthres);
	std::cerr << "done.\n";

	// instantiate bam writer
	::libmaus2::bambam::BamHeader::unique_ptr_type Pbamheader(aligner.getBamHeader(arginfo,std::string(PACKAGE_VERSION)));
	::libmaus2::bambam::BamWriter bamwr(std::cout,*Pbamheader);

	for ( uint64_t irestarg = 1; irestarg < arginfo.restargs.size(); ++irestarg )
	{
		std::string const fqfn = arginfo.getUnparsedRestArg(irestarg);

		std::cerr << "\n" << fqfn << "\n\n";

		// open queries
		libmaus2::aio::InputStreamInstance CIS(fqfn);
		libmaus2::fastx::StreamFastQReaderWrapper queriesIn(CIS);
		libmaus2::fastx::StreamFastQReaderWrapper::pattern_type pattern;

		uint64_t read = 0;
		uint64_t aligned = 0;
		libmaus2::util::Histogram maphist;
		libmaus2::util::Histogram intervalhist;

		// read next query
		while ( queriesIn.getNextPatternUnlocked(pattern) )
		{
			std::string const & t_sid   = pattern.sid;
			std::string const & t_query = pattern.spattern;
			std::string const & t_queryquality = pattern.quality;

			struct ScoreBucketEntry
			{
				uint64_t z;
				libmaus2::alignment::BamLineInfo BLI;

				ScoreBucketEntry() : z(0), BLI()
				{

				}

				ScoreBucketEntry(uint64_t const rz, libmaus2::alignment::BamLineInfo const & rBLI) : z(rz), BLI(rBLI)
				{

				}
			};

			struct ScoreBucketZComparator
			{
				bool operator()(ScoreBucketEntry const & A, ScoreBucketEntry const & B) const
				{
					return A.z < B.z;
				}
			};

			struct ScoreBucketCoordinateComparator
			{
				bool operator()(ScoreBucketEntry const & A, ScoreBucketEntry const & B) const
				{
					if ( A.BLI.refid != B.BLI.refid )
						return A.BLI.refid < B.BLI.refid;
					else if ( A.BLI.pos != B.BLI.pos )
						return A.BLI.pos < B.BLI.pos;
					else
						return false;
				}

				bool operator()(std::vector<ScoreBucketEntry> const & A, std::vector<ScoreBucketEntry> const & B) const
				{
					if ( ! A.size() )
						return true;
					else if ( ! B.size() )
						return false;
					else
						return operator()(A.front(),B.front());
				}
			};

			// alignment vector lock
			libmaus2::parallel::StdMutex BLIlock;
			// all alignments vector
			std::vector<ScoreBucketEntry> allBLIs;

			// maximum fragment start position considered
			int64_t const looplimit = static_cast<int64_t>(t_query.size()) - static_cast<int64_t>(fraglen);

			#if defined(_OPENMP)
			#pragma omp parallel for
			#endif
			for ( int64_t zz = 0; zz <= looplimit; ++zz )
			{
				std::ostringstream sidostr;
				sidostr << t_sid << "_" << read << "_" << std::setw(6) << std::setfill('0') << zz;
				// query sequence id
				std::string const sid = sidostr.str();
				// query string
				std::string const query = t_query.substr(zz,fraglen);
				// quality
				std::string const queryquality = t_queryquality.substr(zz,fraglen);
				// space for mapped query
				libmaus2::autoarray::AutoArray<char> mapped;

				// compute alignments
				std::vector<libmaus2::alignment::BamLineInfo> BLIs;
				aligner.align(sid,query,queryquality,mapped,BLIs);

				// store alignments
				{
					libmaus2::parallel::StdMutex::scope_lock_type SL(BLIlock);

					for ( uint64_t i = 0; i < BLIs.size(); ++i )
						if (
							BLIs[i].getFrontSoftClipping() <= clipthres
							&&
							BLIs[i].getBackSoftClipping() <= clipthres
							&&
							BLIs[i].getErrorRate() <= errrate
						)
							allBLIs.push_back(ScoreBucketEntry(zz,BLIs[i]));
				}
			}

			// sort alignments by z (fragment start)
			std::sort(allBLIs.begin(),allBLIs.end(),ScoreBucketZComparator());
			std::vector< std::vector<ScoreBucketEntry> > scorebuckets;

			// sort alignments into buckets
			uint64_t blow = 0;
			while ( blow != allBLIs.size() )
			{
				// get all alignments for one z value
				uint64_t bhigh = blow;
				while ( bhigh < allBLIs.size() && allBLIs[blow].z == allBLIs[bhigh].z )
					++bhigh;

				// copy alignments to vector BLIs
				uint64_t const zz = allBLIs[blow].z;
				std::vector<libmaus2::alignment::BamLineInfo> BLIs;
				for ( uint64_t i = blow; i < bhigh; ++i )
					BLIs.push_back(allBLIs[i].BLI);

				blow = bhigh;

				// sort alignments into buckets, bucket is chosen by considering
				// z value and mapping coordinate
				for ( uint64_t i = 0; i < BLIs.size(); ++i )
				{
					int64_t ins = -1;
					for ( uint64_t j = 0; j < scorebuckets.size(); ++j )
						if (
							BLIs[i].refid == scorebuckets[j].back().BLI.refid
							&&
							absdif(
								// difference in mapping position
								absdif(BLIs[i].pos,scorebuckets[j].back().BLI.pos)
								,
								// difference in position on read
								absdif(zz,scorebuckets[j].back().z)
							)
							<= difthres
						)
							ins = j;

					if ( ins >= 0 )
						scorebuckets[ins].push_back(ScoreBucketEntry(zz,BLIs[i]));
					else
						scorebuckets.push_back(std::vector<ScoreBucketEntry>(1,ScoreBucketEntry(zz,BLIs[i])));
				}

				#if 0
				if ( zz % 512 == 0 )
					std::cerr << "zz=" << zz << "/" << (t_query.size()-fraglen+1) << " number of buckets " << scorebuckets.size() << std::endl;
				#endif

				// bam encoding for fragment
				if ( BLIs.size() )
				{
					aligned++;

					if ( BLIs.size() > 1 )
					{
						BLIs[0].encode(bamwr);
						for ( uint64_t i = 1; i < BLIs.size(); ++i )
						{
							BLIs[i].flags |= libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY;
							BLIs[i].encode(bamwr);
						}
					}
					else if ( BLIs.size() == 1 )
					{
						BLIs[0].encode(bamwr);
					}
				}
				else
				{
					#if 0
					std::cerr << "[V] no alignments for " << sid << " of length " << query.size() << std::endl;
					#endif
				}

				// update mapping count histogram
				maphist(BLIs.size());
			}

			// reverse order if we mapped in the opposite direction
			for ( uint64_t i = 0; i < scorebuckets.size(); ++i )
				if (
					scorebuckets[i].size() &&
					scorebuckets[i].front().BLI.pos > scorebuckets[i].back().BLI.pos
				)
					std::reverse(scorebuckets[i].begin(),scorebuckets[i].end());

			// sort by mapping coordinate
			std::sort(scorebuckets.begin(),scorebuckets.end(),ScoreBucketCoordinateComparator());

			typedef std::pair<uint64_t,uint64_t> upair;
			std::vector<upair> I;
			uint64_t ilow = 0;

			std::cerr << "\n" << std::string(1,' ') << t_sid << " (" << t_query.size() << ")" << " score buckets " << scorebuckets.size() << std::endl;

			// sort alignments into region buckets
			uint64_t scorebucketdist = 10000;
			while ( ilow != scorebuckets.size() )
			{
				uint64_t ihigh = ilow+1;
				while (
					(ihigh < scorebuckets.size()) &&
					// same ref seq
					(scorebuckets[ihigh-1].back().BLI.refid == scorebuckets[ihigh].front().BLI.refid) &&
					// no more than scorebucketdist distance in raw mapping position
					(
						absdif(scorebuckets[ihigh-1].back().BLI.pos,scorebuckets[ihigh].front().BLI.pos)
						<= scorebucketdist
					)
				)
					++ihigh;

				#if 0
				if (
					Pbamheader->getRefIDName(scorebuckets[ilow].front().BLI.refid) == std::string("chr1") &&
					scorebuckets[ilow].front().BLI.pos >= 155000000ll &&
					scorebuckets[ihigh-1].back().BLI.pos <= 170000000ll
				)
				#endif
				{
					I.push_back(upair(ilow,ihigh));
				}

				ilow = ihigh;
			}

			//plot "ABC10_2_1_000044079600_C10_newbler_large_blastn_cpatch_gap5_consensus.fq_merged_4_0.txt" with lines title "1",
			//     "ABC10_2_1_000044079600_C10_newbler_large_blastn_cpatch_gap5_consensus.fq_merged_4_1.txt" with lines title "2",
			//     "ABC10_2_1_000044079600_C10_newbler_large_blastn_cpatch_gap5_consensus.fq_merged_4_2.txt" with lines title "3"

			if ( I.size() )
			{
				std::vector<ZSortEntry> ZSEV;

				std::ostringstream plotstr;
				plotstr << "set terminal postscript eps color\n";
				plotstr << "set style line 1 lt 1 lw 2 pt 1 ps 1 linecolor rgb \"black\"\n";
				plotstr << "set style line 2 lt 2 lw 2 pt 2 ps 1 linecolor rgb \"red\"\n";
				plotstr << "set style line 3 lt 3 lw 2 pt 3 ps 1 linecolor rgb \"green\"\n";
				plotstr << "set style line 4 lt 4 lw 2 pt 4 ps 1 linecolor rgb \"blue\"\n";
				plotstr << "plot ";

				double minscore = std::numeric_limits<double>::max();
				double maxscore = -std::numeric_limits<double>::max();
				double zlow = std::numeric_limits<double>::max();
				double zhigh = -std::numeric_limits<double>::max();
				for ( uint64_t j = 0; j < I.size(); ++j )
					for ( uint64_t i = I[j].first; i < I[j].second; ++i )
						for ( uint64_t k = 0; k < scorebuckets[i].size(); ++k )
						{
							double const score = scorebuckets[i][k].BLI.score;
							minscore = std::min(minscore,score);
							maxscore = std::max(maxscore,score);
							double const z = scorebuckets[i][k].z;
							zlow = std::min(zlow,z);
							zhigh = std::max(zhigh,z);
						}

				plotstr << "[" << zlow << " to " << zhigh << "] " << "[" << minscore << " to " << maxscore+30 << "]";
				uint64_t p = 0;
				uint64_t pfnid = 0;

				std::vector<std::string> datafns;
				for ( uint64_t j = 0; j < I.size(); ++j )
				{
					std::map<int64_t,double> scoremap;

					std::cerr << std::string(2,' ') << std::string(80,'-') << std::endl;

					for ( uint64_t i = I[j].first; i < I[j].second; ++i )
					{
						double scoresum = 0;
						for ( uint64_t j = 0; j < scorebuckets[i].size(); ++j )
							scoresum += static_cast<double>(scorebuckets[i][j].BLI.score) / scorebuckets[i][j].BLI.seq.size();

						uint64_t offthres = fraglen;
						uint64_t klow = 0;
						while ( klow != scorebuckets[i].size() )
						{
							uint64_t khigh = klow+1;
							while ( khigh != scorebuckets[i].size() && absdif(scorebuckets[i][khigh].z,scorebuckets[i][khigh-1].z) <= offthres )
								++khigh;

							std::cerr << "klow=" << klow << " khigh=" << khigh << std::endl;

							if ( klow != 0 )
								std::cerr << scorebuckets[i][klow-1].z << "," << scorebuckets[i][klow].z << std::endl;

							uint64_t const zlow  = std::min(scorebuckets[i][klow].z,scorebuckets[i][khigh-1].z);
							uint64_t const zhigh = std::max(scorebuckets[i][klow].z,scorebuckets[i][khigh-1].z);
							uint64_t const zend = zhigh + fraglen;
							bool const zreverse = scorebuckets[i][klow].z > scorebuckets[i][khigh-1].z;
							std::string const refq = t_query.substr(zlow,zend-zlow);

							uint64_t const seq    = scorebuckets[i][klow].BLI.refid;
							uint64_t const seqpos = scorebuckets[i][klow].BLI.pos;
							uint64_t const seqlen = scorebuckets[i][khigh-1].BLI.pos + fraglen - seqpos;
							std::string refsub = aligner.index.getTextUnmapped(aligner.ISA[aligner.seqstarts[2*seq] + seqpos],seqlen);
							if ( zreverse )
								refsub = libmaus2::fastx::reverseComplementUnmapped(refsub);

							libmaus2::lcs::MetaLocalEditDistance< ::libmaus2::lcs::diag_del_ins > LED;
							libmaus2::lcs::LocalEditDistanceResult LEDR = LED.process(refsub.begin(),refsub.size(),refq.begin(),refq.size(),
								::std::floor((errrate * std::max(refsub.size(),refq.size()))+0.5)
							);

							if ( LEDR.getErrorRate() <= errrate )
							{
								libmaus2::lcs::LocalAlignmentPrint::printAlignmentLines(std::cerr,refsub,refq,80,LED.ta,LED.te,LEDR);

								std::cerr << "[V] error rate " << LEDR.getErrorRate()
									<< " zlow=" << zlow << " zend=" << zend << " zreverse=" << zreverse
									<< " seq=" << Pbamheader->getRefIDName(seq) << " seqpos=" << seqpos << " seqlen=" << seqlen
									<< std::endl;

								std::vector<libmaus2::bambam::BamFlagBase::bam_cigar_ops> cigopvec;
								libmaus2::lcs::LocalAlignmentTraceContainer const & trace = LED.getTrace();
								for (
									libmaus2::lcs::LocalAlignmentTraceContainer::step_type const * tc = trace.ta;
									tc != trace.te;
									++tc )
								{
									switch ( *tc )
									{
										case libmaus2::lcs::LocalAlignmentTraceContainer::STEP_MATCH:
											cigopvec.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL);
											break;
										case libmaus2::lcs::LocalAlignmentTraceContainer::STEP_MISMATCH:
											cigopvec.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF);
											break;
										case libmaus2::lcs::LocalAlignmentTraceContainer::STEP_INS:
											cigopvec.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS);
											break;
										case libmaus2::lcs::LocalAlignmentTraceContainer::STEP_DEL:
											cigopvec.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL);
											break;
										default:
											break;
									}
								}

								uint64_t seqclipleft = zreverse ? LEDR.a_clip_right : LEDR.a_clip_left;
								uint64_t seqclipright = zreverse ? LEDR.a_clip_left : LEDR.a_clip_right;

								ZSEV.push_back(
									ZSortEntry(
										zlow+LEDR.b_clip_left,
										zend-LEDR.b_clip_right,
										zreverse,
										seq,
										seqpos+seqclipleft,seqpos+seqlen-seqclipright,
										cigopvec
									)
								);

								std::ostringstream pfnstr;
								pfnstr << fqfn << "_" << t_sid << "_" << j << "_" << (pfnid++) << "_" << (zend-zlow)-(LEDR.b_clip_left+LEDR.b_clip_right) << ".pairwise.fasta";
								std::string const pfn = pfnstr.str();

								libmaus2::aio::OutputStreamInstance PFNCOS(pfn);
								PFNCOS << ">" << Pbamheader->getRefIDName(seq) << "_" << seqpos+seqclipleft << "_" << seqpos+seqlen-seqclipright << "_" << zreverse << '\n';
								PFNCOS << refsub.substr(LEDR.a_clip_left,refsub.size()-(LEDR.a_clip_left+LEDR.a_clip_right)) << '\n';
								PFNCOS << ">" << t_sid << "_" << zlow+LEDR.b_clip_left << "_" << zend-LEDR.b_clip_right << '\n';
								PFNCOS << refq.substr(LEDR.b_clip_left,refq.size()-(LEDR.b_clip_left+LEDR.b_clip_right)) << '\n';
								PFNCOS.flush();
								//PFNCOS.close();
							}

							klow = khigh;
						}

						std::cerr << std::string(3,' ') << "[" << i << "]="
							<< scorebuckets[i].size()
							<< "\t"
							// could be descending
							<< scorebuckets[i].front().z << "," << scorebuckets[i].back().z
							<< "\t"
							// ascending coordinates as sorted that way
							<< "("
							<< Pbamheader->getRefIDName(scorebuckets[i].front().BLI.refid)
							<< "," << scorebuckets[i].front().BLI.pos
							<< "," << scorebuckets[i].back().BLI.pos
							<< ")"
							<< "\t"
							<< scoresum
							<< std::endl;

						std::vector<ScoreBucketEntry> const & V = scorebuckets[i];
						for ( uint64_t k = 0; k < V.size(); ++k )
							scoremap[ V[k].z ] = V[k].BLI.score;
					}

					typedef std::pair<int64_t,double> idtype;
					std::vector<idtype> const ID(scoremap.begin(),scoremap.end());
					std::vector<upair> IDintervals;

					uint64_t idlow = 0;
					while ( idlow < ID.size() )
					{
						uint64_t idhigh = idlow+1;
						while ( idhigh < ID.size() && ID[idhigh-1].first+1 == ID[idhigh].first )
							++idhigh;

						IDintervals.push_back(upair(idlow,idhigh));

						idlow = idhigh;
					}

					for ( uint64_t idi = 0; idi < IDintervals.size(); ++idi )
					{
						std::ostringstream gplfnstr;
						gplfnstr << fqfn << "_" << t_sid << "_" << j << "_" << idi << ".txt";
						std::string const gplfn = gplfnstr.str();
						datafns.push_back(gplfn);

						plotstr << ((p++)==0?" ":",") << "\"" << gplfn << "\" with lines title \"";

						if ( idi == 0 )
							plotstr
								<< Pbamheader->getRefIDName(scorebuckets[I[j].first].front().BLI.refid)
								<< "["
								<< scorebuckets[I[j].first].front().BLI.pos << "-"
								<<
								scorebuckets[I[j].second-1].back().BLI.pos+fraglen
								<< "]";

						plotstr
							<< "\" ls " << (j+1);

						libmaus2::aio::OutputStreamInstance gplCOS(gplfn);
						for ( uint64_t idii = IDintervals[idi].first; idii < IDintervals[idi].second; ++idii )
						{
							gplCOS << ID[idii].first << "\t" << ID[idii].second << "\n";
						}
						gplCOS.flush();
						//gplCOS.close();
					}
				}

				plotstr << "\n";
				std::string const plotseq = plotstr.str();

				std::ostringstream gplfnstr;
				gplfnstr << fqfn << "_" << t_sid << ".gpl";
				std::string const gplfn = gplfnstr.str();
				std::ostringstream epsfnstr;
				epsfnstr << fqfn << "_" << t_sid << ".eps";
				std::string const epsfn = epsfnstr.str();

				{
				libmaus2::aio::OutputStreamInstance COS(gplfn);
				COS << plotseq;
				COS.flush();
				}
				// COS.close();

				std::ostringstream plotcmdstr;
				plotcmdstr << "gnuplot < " << gplfn << " >" << epsfn;
				std::string const plotcmd = plotcmdstr.str();
				int r = system(plotcmd.c_str());

				if ( r != EXIT_SUCCESS )
				{
					std::cerr << plotcmd << " failed" << std::endl;
				}

				std::ostringstream epstopdfstr;
				epstopdfstr << "epstopdf " << epsfn;
				std::string const epstopdf = epstopdfstr.str();
				r = system(epstopdf.c_str());

				if ( r != EXIT_SUCCESS )
				{
					std::cerr << epstopdf << " failed" << std::endl;
				}

				std::sort(ZSEV.begin(),ZSEV.end());

				uint64_t zzid = 0;

				uint64_t zzlow = 0;
				while ( zzlow != ZSEV.size() )
				{
					uint64_t zzhigh = zzlow;
					while ( zzhigh != ZSEV.size() && ZSEV[zzhigh].zreverse == ZSEV[zzlow].zreverse )
						++zzhigh;

					if ( zzlow == 0 && zzhigh != ZSEV.size() )
						std::cerr << "[V] warning, contig maps in both directions." << std::endl;

					std::deque<ZSortEntry> subvec(ZSEV.begin()+zzlow,ZSEV.begin()+zzhigh);

					std::cerr << std::string(80,'+') << std::endl;
					for ( uint64_t i = 0; i < subvec.size(); ++i )
						std::cerr << subvec[i] << std::endl;

					while ( subvec.size() )
					{
						// collect intervals starting at the same z value
						uint64_t iend = 0;
						uint64_t minend = subvec[0].getZEnd();
						while (
							iend != subvec.size() &&
							subvec[iend].getZStart() == subvec[0].getZStart()
						)
							++iend;

						// shrink to overlap
						for ( uint64_t i = iend; i < subvec.size(); ++i )
							if ( subvec[i].getZStart() >= minend )
							{
								// std::cerr << "Break on zstart=" << subvec[i].getZStart() << std::endl;
								break;
							}
							else
							{
								minend = std::min(minend,subvec[i].getZStart());
								// std::cerr << "Set minend to " << minend << " for " << subvec[i] << std::endl;
							}

						// std::cerr << "minend=" << minend << std::endl;

						std::vector < ZSortEntry > writevec;
						std::cerr << std::string(80,'-') << std::endl;
						for ( uint64_t i = 0; i < iend; ++i )
						{
							std::pair < ZSortEntry, ZSortEntry > ZP = subvec[i].zsplit(minend);

							// std::cerr << ZP.first << std::endl;

							uint64_t const zlow  = ZP.first.zstart;
							uint64_t const zend = ZP.first.zend;
							bool const zreverse = ZP.first.zreverse;
							std::string const refq = t_query.substr(zlow,zend-zlow);

							uint64_t const seq    = ZP.first.seq;
							uint64_t const seqpos = ZP.first.seqstart;
							uint64_t const seqlen = ZP.first.seqend-ZP.first.seqstart;
							std::string refsub = aligner.index.getTextUnmapped(aligner.ISA[aligner.seqstarts[2*seq] + seqpos],seqlen);
							if ( zreverse )
								refsub = libmaus2::fastx::reverseComplementUnmapped(refsub);

							libmaus2::lcs::MetaLocalEditDistance< ::libmaus2::lcs::diag_del_ins > LED;
							libmaus2::lcs::LocalEditDistanceResult LEDR = LED.process(refsub.begin(),refsub.size(),refq.begin(),refq.size(),
								::std::floor((errrate * std::max(refsub.size(),refq.size()))+0.5)
							);
							// libmaus2::lcs::LocalAlignmentPrint::printAlignmentLines(std::cerr,refsub,refq,80,LED.ta,LED.te,LEDR);

							if ( zreverse )
							{
								ZP.first.seqend -= LEDR.a_clip_left;
								ZP.first.seqstart += LEDR.a_clip_right;
							}
							else
							{
								ZP.first.seqstart += LEDR.a_clip_left;
								ZP.first.seqend -= LEDR.a_clip_right;
							}

							writevec.push_back(ZP.first);

							subvec[i] = ZP.second;
						}

						std::ostringstream mafnstr;
						mafnstr << fqfn << "_" << t_sid
							<< "_" << std::setw(6) << std::setfill('0') << (zzid++) << std::setw(0)
							<< "_" << std::setw(6) << std::setfill('0') << writevec[0].zstart << std::setw(0)
							<< "_" << std::setw(6) << std::setfill('0') << writevec[0].zend << std::setw(0)
							<< ".multi.fasta";
						libmaus2::aio::OutputStreamInstance COS(mafnstr.str());

						{
							uint64_t const zlow = writevec[0].zstart;
							uint64_t const zend = writevec[0].zend;
							std::string const refq = t_query.substr(zlow,zend-zlow);
							COS << ">contig_" << zlow << "_" << zend << '\n';
							COS << refq << '\n';
						}


						for ( uint64_t i = 0; i < writevec.size(); ++i )
						{
							ZSortEntry const & ZP = writevec[i];

							std::cerr << "[" << i << "]=" << writevec[i] << std::endl;

							bool const zreverse = ZP.zreverse;
							uint64_t const seq    = ZP.seq;
							uint64_t const seqpos = ZP.seqstart;
							uint64_t const seqlen = ZP.seqend-ZP.seqstart;
							std::string refsub = aligner.index.getTextUnmapped(aligner.ISA[aligner.seqstarts[2*seq] + seqpos],seqlen);
							if ( zreverse )
								refsub = libmaus2::fastx::reverseComplementUnmapped(refsub);

							COS << ">region_" << Pbamheader->getRefIDName(seq) << "_" << seqpos << "_" << seqpos+seqlen << '\n';
							COS << refsub << '\n';
						}

						COS.flush();
						// COS.close();

						while ( subvec.size() && (subvec.front().zend==subvec.front().zstart) )
							subvec.pop_front();

						std::sort(subvec.begin(),subvec.end());
					}

					zzlow = zzhigh;
				}
			}

			read += 1;
		}

		// std::cerr << read << "\t" << aligned << "\t" << static_cast<double>(aligned)/read << std::endl;

		{
			std::cerr << "\n map histogram" << std::endl;
			std::ostringstream ostr;
			maphist.print(ostr);
			std::stringstream istr(ostr.str());
			while ( istr )
			{
				std::string line;
				std::getline(istr,line);
				if ( line.size() )
					std::cerr << " " << line << std::endl;
			}
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc, argv);

		if ( arginfo.restargs.size() < 1 )
		{
			std::cerr << "[E] usage: " << argv[0] << " <indexprefix>" << std::endl;
			return EXIT_FAILURE;
		}

		std::string const mode = arginfo.getValue<std::string>("mode","hwt");

		if ( mode == "hwt" )
		{
			return alignalt<libmaus2::lf::ImpCompactHuffmanWaveletLF>(arginfo,".hwt");
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unknown mode " << mode << std::endl;
			lme.finish();
			throw lme;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
