### Build (old way, unstable)
I don't use this anymore as it was unreliable and not forward compatible.

```sh
./prepare_conf.sh

./configure --enable-maintainer-mode --enable-debug

make

make check

make install
```

may give some errors on the newest version of boost.

