#! /usr/bin/env python
# encoding: utf-8

VERSION = '0.0.1'
APPNAME = 'ripples'


def options(opt):
    opt.load('compiler_cxx')
    opt.load('waf_unit_test')
    opt.load('sphinx', tooldir='waftools')

    cfg_options = opt.get_option_group('Configuration options')

    opt.load('trng4', tooldir='waftools')
    opt.load('json', tooldir='waftools')
    opt.load('spdlog', tooldir='waftools')

    cfg_options.add_option(
        '--openmp-root', action='store', default='/usr',
        help='root directory of the installation of openmp')

    cfg_options.add_option(
        '--enable-mpi', action='store_true', default=False,
        help='enable openmpi implementation')

    opt.load('mpi', tooldir='waftools')


def configure(conf):
    conf.load('compiler_cxx')
    conf.load('waf_unit_test')
    conf.load('sphinx', tooldir='waftools')

    conf.env.CXXFLAGS += ['-std=c++17', '-O3', '-march=native', '-pipe']

    conf.load('spdlog', tooldir='waftools')
    conf.load('json', tooldir='waftools')
    conf.load('trng4', tooldir='waftools')

    conf.check_cxx(cxxflags=['-fopenmp'], ldflags=['-fopenmp'],
                   libpath=['{0}/lib/'.format(conf.options.openmp_root)],
                   uselib_store='OpenMP')

    if conf.options.enable_mpi:
        conf.load('mpi', tooldir='waftools')


def build(bld):
    directories = ['3rd-party', 'include', 'tools', 'test']

    bld.recurse(directories)


def build_docs(bld):
    if bld.env.ENABLE_DOCS:
        bld(features='sphinx', sources='docs')
    else:
        bld.fatal('Please configure with --enable-docs')

from waflib import Build
class docs(Build.BuildContext):
    fun = 'build_docs'
    cmd = 'docs'
