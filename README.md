alignalt
========

Program aligning contigs to reference regions for detecting matching stretches

Compilation of alignalt
---------------------------

alignalt needs libmaus2 [https://github.com/gt1/libmaus2].
When libmaus2 is installed in ${LIBMAUS2PREFIX}
then alignalt can be compiled and installed in ${HOME}/alignalt using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUS2PREFIX} \
		--prefix=${HOME}/alignalt
	- make install
