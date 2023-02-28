Describe how to configure the workflow (using config.yaml and maybe additional files).
All of them need to be present with example entries inside of the config folder.

you need to install the latest version of foldseek
the binaries on conda are out of date so its not included in the snakemake env

# static Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
# static Linux SSE4.1 build (check using: cat /proc/cpuinfo | grep sse4_1)
wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
# static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

One of these should do it. foldseek path is set in the workflow at the top. it defaults to project_dir/foldseek/bin/foldseek.
