import subprocess
import sys

importfromfile = __import__(sys.argv[1],  fromlist=['boost_versions',  'nvcc_versions',  'gcc_versions',  'cuda_arch'])
boost_versions = getattr(importfromfile,  'boost_versions')
nvcc_versions = getattr(importfromfile,  'nvcc_versions')
gcc_versions = getattr(importfromfile,  'gcc_versions')
cuda_arch = getattr(importfromfile,  'cuda_arch')

applications = ['LandauGaugeDP',  'CoulombGaugeDP',  'LandauGaugeSP',  'CoulombGaugeSP']
cpp11ornot = { 'c++11': '-DCULGT_ALLOW_C++11=ON',  'NOc++11': '-DCULGT_ALLOW_C++11=OFF'}

otheroptions = '-DDISABLE_AUTOTUNE=ON'

subprocess.Popen(['mkdir',  'build']).wait()

for boost_version in boost_versions.keys():
    for nvcc_version in nvcc_versions.keys():
        for gcc_version in gcc_versions.keys():
            for cpp11 in cpp11ornot.keys():
                subprocess.Popen(['make',  'clean'],  cwd='build' ).wait()
                subprocess.Popen(['rm',  'CMakeCache.txt'],  cwd='build' ).wait()
                subprocess.Popen(['cmake',  '../../../lib/gaugefixing',   cpp11ornot[cpp11], otheroptions,  '-DCUDA_ARCHITECTURE='+cuda_arch,'-DBOOST_ROOT='+boost_versions[boost_version], '-DCUDA_HOST_COMPILER='+gcc_versions[gcc_version], '-DCUDA_TOOLKIT_ROOT_DIR='+nvcc_versions[nvcc_version] ],  cwd='build' ).wait()
                for app in applications:
                    subprocess.Popen(['make', app],  cwd='build' ).wait()
                    subprocess.Popen(['mv', app,  app+'_' + nvcc_version + '_' + gcc_version + '_' + boost_version + '_' + cpp11],  cwd='build' ).wait()

