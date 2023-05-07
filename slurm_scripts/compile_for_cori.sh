pipenv shell

module unload PrgEnv-intel
module load PrgEnv-gnu

module use /global/common/software/m3169/cori/modulefiles
module load openmpi

conan create conan/waf-generator user/stable
conan create conan/trng user/stable
conan create conan/metall user/stable
conan create conan/memkind user/stable
conan install --install-folder build .
