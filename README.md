alignalt
========

Program aligning contigs to reference regions for detecting matching stretches

Compilation of alignalt
---------------------------

alignalt needs libmaus [https://github.com/gt1/libmaus].
When libmaus is installed in ${LIBMAUSPREFIX}
then alignalt can be compiled and installed in ${HOME}/alignalt using

	- autoreconf -i -f
	- ./configure --with-libmaus=${LIBMAUSPREFIX} \
		--prefix=${HOME}/alignalt
	- make install
