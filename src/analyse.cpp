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
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/fastx/FastAStreamSet.hpp>

struct AlignmentComparator
{
	bool operator()(libmaus2::bambam::BamAlignment::shared_ptr_type const & A, libmaus2::bambam::BamAlignment::shared_ptr_type const & B)
	{
		std::string const name_a = A->getName();
		std::string const name_b = B->getName();

		if ( name_a != name_b )
			return name_a < name_b;

		int32_t const ref_len_a = A->isMapped() ? A->getReferenceLength() : -1;
		int32_t const ref_len_b = B->isMapped() ? B->getReferenceLength() : -1;

		// descending
		if ( ref_len_a != ref_len_b )
			return ref_len_a > ref_len_b;

		int32_t const ref_id_a = A->getRefID();
		int32_t const ref_id_b = B->getRefID();

		if ( ref_id_a != ref_id_b )
			return ref_id_a < ref_id_b;

		int const a_is_rev = A->isReverse();
		int const b_is_rev = B->isReverse();

		if ( a_is_rev != b_is_rev )
			return a_is_rev < b_is_rev;

		return false;
	}
};

struct PosComparator
{
	bool operator()(libmaus2::bambam::BamAlignment::shared_ptr_type const & A, libmaus2::bambam::BamAlignment::shared_ptr_type const & B)
	{
		return A->getPos() < B->getPos();
	}
};

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		if ( ! arginfo.hasArg("reference") )
		{
			std::cerr << "[E] missing reference key" << std::endl;
			return EXIT_FAILURE;
		}

		std::string const reffn = arginfo.getUnparsedValue("reference",std::string());
		std::string const seqid = arginfo.getUnparsedValue("seqid",std::string());

		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type Pdec(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
		libmaus2::bambam::BamAlignmentDecoder & decoder = Pdec->getDecoder();
		libmaus2::bambam::BamHeader const & header = decoder.getHeader();

		assert (header.getNumRef()) ;
		std::string const ref = libmaus2::fastx::FastAStreamSet::getStreamAsString(reffn,header.getRefIDName(0));
		// std::cerr << ref << std::endl;

		//std::cerr << header.text << std::endl;
		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type > Valgn;
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigopin;
		libmaus2::autoarray::AutoArray<char> readdata;
		libmaus2::bambam::BamAlignment::D_array_type T;
		libmaus2::bambam::BamFormatAuxiliary formataux;

		while ( decoder.readAlignment() )
		{
			libmaus2::bambam::BamAlignment const & algn = decoder.getAlignment();
			if ( algn.isMapped() )
			{
				assert ( algn.getRefID() == 0 );

				libmaus2::bambam::BamAlignment::shared_ptr_type salgn(algn.sclone());
				uint64_t const numcig = libmaus2::bambam::BamAlignmentDecoderBase::recalculateCigar(
					salgn->D.begin(),
					ref.begin() + salgn->getPos(),
					cigopin,
					readdata
				);
				salgn->replaceCigarString(cigopin,numcig,T);

				#if 0
				salgn->formatAlignment(std::cerr,header,formataux);
				std::cerr << std::endl;
				#endif

				Valgn.push_back(salgn);
			}
		}

		std::sort(Valgn.begin(),Valgn.end(),AlignmentComparator());

		uint64_t low = 0;
		uint64_t out = 0;
		while ( low < Valgn.size() )
		{
			uint64_t high = low+1;
			while (
				high < Valgn.size() &&
				strcmp(Valgn[low]->getName(),Valgn[high]->getName()) == 0
			)
			{
				++high;
			}

			Valgn[out++] = Valgn[low];

			low = high;
		}

		Valgn.resize(out);

		for ( uint64_t i = 0; i < Valgn.size(); ++i )
		{
			libmaus2::bambam::BamAlignment const & algn = *Valgn[i];
			uint64_t const po = algn.getPos();
			uint64_t const rl = algn.getReferenceLength();
		}

		out = 0;
		for ( uint64_t i = 0; i < Valgn.size(); ++i )
		{
			bool overlapfree = true;
			for ( uint64_t j = 0; j < out; ++j )
			{
				int64_t const prev_start = Valgn[j]->getPos();
				int64_t const prev_end   = prev_start + static_cast<int64_t>(Valgn[j]->getReferenceLength()) - 1;
				int64_t const n_start = Valgn[i]->getPos();
				int64_t const n_end   = n_start + static_cast<int64_t>(Valgn[i]->getReferenceLength()) - 1;
				libmaus2::math::IntegerInterval<int64_t> prev(prev_start,prev_end);
				libmaus2::math::IntegerInterval<int64_t> nintv(n_start,n_end);
				libmaus2::math::IntegerInterval<int64_t> ints = libmaus2::math::IntegerInterval<int64_t>::intersection(prev,nintv);

				// std::cerr << prev << " " << nintv << " " << ints << std::endl;

				if ( ! ints.isEmpty() )
					overlapfree = false;
			}

			if ( overlapfree )
				Valgn[out++] = Valgn[i];
		}

		Valgn.resize(out);

		std::sort(Valgn.begin(),Valgn.end(),PosComparator());

		libmaus2::autoarray::AutoArray<uint64_t> cigstats(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CTHRES);
		std::ostringstream intvstr;
		std::ostringstream datastr;
		uint64_t prevend = 0;
		for ( uint64_t i = 0; i < Valgn.size(); ++i )
		{
			libmaus2::bambam::BamAlignment const & algn = *Valgn[i];
			uint64_t const po = algn.getPos();
			uint64_t const rl = algn.getReferenceLength();
			libmaus2::math::IntegerInterval<int64_t> nintv(po,po+rl-1);

			std::string data = algn.getRead();

			std::string const frontclipped = data.substr(0,algn.getFrontSoftClipping());
			std::string const backclipped = data.substr(data.size() - algn.getBackSoftClipping());

			data = data.substr(algn.getFrontSoftClipping());
			data = data.substr(0,data.size() - algn.getBackSoftClipping());

			intvstr << algn.getName() << ":" << data << ":" << frontclipped << ":" << backclipped << ":" << nintv << ";";

			algn.getCigarStats(cigstats,false);

			if ( po != prevend )
			{
				datastr << ref.substr(prevend,po-prevend);
			}

			datastr << data;

			prevend = po + rl;
		}

		if ( prevend != ref.size() )
			datastr << ref.substr(prevend);

		double const ident = libmaus2::bambam::BamAlignmentDecoderBase::getIdentityFractionFromCigarStats(cigstats);
		std::cerr << seqid << "\t" << reffn << "\t" << ident << "\t" << cigstats[libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL] << "\t"
			<< cigstats[libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF] << "\t"
			<< intvstr.str()
			<< std::endl;

		std::cout << ">" << seqid << "\n" << datastr.str() << std::endl;

		#if 0
		for ( uint64_t i = 0; i < cigstats.size(); ++i )
			std::cerr << i << " " << cigstats[i] << std::endl;
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
