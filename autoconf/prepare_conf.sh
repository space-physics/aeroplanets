#!/usr/bin/env bash
aclocal --force
autoheader
libtoolize --force
#touch NEWS README AUTHORS ChangeLog
automake --add-missing --copy --force-missing --gnu
autoconf

autoreconf -isvf

