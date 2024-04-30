GreediRIS
*******

This repository contains the software application described in *Scalable Influence Maximization using Distributed Streaming Maximum Cover*.

Quickstart with Conan
=====================

First of all we need to set up the Python environment needed.

.. code-block:: shell

   $ pip install --user pipenv
   $ pipenv --three
   $ pipenv install
   $ pipenv shell

Then we need to install dependencies:

.. code-block:: shell

   $ conan create conan/waf-generator user/stable
   $ conan create conan/trng user/stable
   $ conan create conan/metall user/stable
   $ conan create conan/memkind user/stable
   $ conan install --install-folder build .

To enable Memkind or Metall please replace the conan install command with one of:

.. code-block:: shell

   $ conan install --install-folder build . -o memkind=True
   $ conan install --install-folder build . -o metall=True


Now we are ready to configure and build GreediRIS:

.. code-block:: shell

   $ ./waf configure --enable-mpi build_release


Build Instructions
==================

This project uses `WAF <https://waf.io>`_ as its build system.  Building GreediRIS
is a two-step process: configure the project and build the tools.  Before
attempting to build, be sure to have the following dependencies installed:

- A compiler with C++14 support and OpenMP support.
- `Spdlog <https://github.com/gabime/spdlog>`_
- `JSON <https://github.com/nlohmann/json>`_
- `TRNG4 <https://github.com/rabauke/trng4>`_
- An MPI library

The configure step can be invoked with:

.. code-block:: shell

   $ ./waf configure --enable-mpi build_release

The build system offers options that can be used to help the configuration step
locate dependencies (e.g., they are installed in unconventional paths).  A
complete list of the options can be obtained with:

.. code-block:: shell

   $ ./waf configure --help

For more detailed instruction, please read :ref:`build:Step By Step Build
Instructions`.

The tools compiled can be found under ``build/release/tools/``.  A complete set of
command line options can be obtained through:

.. code-block:: shell

   $ ./build/release/tools/<tool_name> --help

Running GreediRIS
================

GreediRIS can be run with ``build/release/tools/mpi-greedi-im``. Running ``build/release/tools/mpi-greedi-im -h`` will provide the following information; 

.. code-block::
   
   Usage: ./build/release/tools/mpi-greedi-im [OPTIONS]

   Options:
      -h,--help                   Print this help message and exit
      

   Input Options:
      -i,--input-graph TEXT REQUIRED
                                    The input file with the edge-list.
      --reload-binary             Reload a graph from binary input
      -u,--undirected             The input graph is undirected
      -w,--weighted               The input graph is weighted
      --distribution TEXT         The distribution to be used (uniform|normal) to generate weights
      --mean FLOAT                The mean for the normal distribution
      --variance FLOAT            The variance for the normal distribution
      --scale-factor FLOAT        Scaling Factor for the generated weights
      --disable-renumbering       Load the graph as is from the input.


   Algorithm Options:
      -k,--seed-set-size UINT REQUIRED
                                    The size of the seed set.
      -p,--parallel               Trigger the parallel implementation
      -d,--diffusion-model TEXT REQUIRED
                                    The diffusion model to use (LT|IC)
      -e,--epsilon FLOAT REQUIRED The size of the seed set.


   Streaming-Engine Options:
      --dump-sampling-data BOOLEAN
                                    Output all sampling data to your output file
      --run-streaming BOOLEAN     Run max-k-cover within a streaming algorithm. False by default.
      --epsilon-2 FLOAT           Set the error parameter for the streaming step. Default of 0.13 to acheive approximation garuntee of 21%
      --alpha FLOAT               Set the fraction of local seeds to send to the final selection step, defaults to 1


   Output Options:
      -o,--output TEXT            The file name of the log.


GreediRIS Team
============

- `Reet Barik | WSU <reet.barik@wsu.edu>`_
- `Wade Cappa | WSU <wade.cappa@wsu.edu>`_
- `S M Ferdous | PNNL <sm.ferdous@pnnl.gov>`_
- `Marco Mintutoli | PNNL <marco.minutoli@pnnl.gov>`_
- `Mahantesh Halappanavar | PNNL, WSU <mahantesh.halappanavar@pnnl.gov>`_
- `Ananth Kalyanaraman | WSU, PNNL <ananth@wsu.edu>`_

This software was produced in collaboration between authors from Washington State University Pullman, and Pacific Northwest National Laboratory Richland. 

Disclamer Notice
================

This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or any
information, apparatus, product, software, or process disclosed, or represents
that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

.. raw:: html

   <div align=center>
   <pre style="align-text:center">
   PACIFIC NORTHWEST NATIONAL LABORATORY
   operated by
   BATTELLE
   for the
   UNITED STATES DEPARTMENT OF ENERGY
   under Contract DE-AC05-76RL01830
   </pre>
   </div>
