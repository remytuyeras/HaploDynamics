from distutils.core import setup

setup(
  name = 'HaploDynamics',
  packages = ['HaploDynamics'],
  version = 'v0.1-beta.0',
  license='gpl-3.0',
  description = 'A python library to develop genomic data simulators',
  author = 'Remy Tuyeras',
  author_email = 'rtuyeras@gmail.com',
  url = 'https://github.com/remytuyeras/HaploDynamics',
  download_url = 'https://github.com/remytuyeras/HaploDynamics/archive/refs/tags/v0.1-beta.0.tar.gz',
  keywords = ['Simulator', 'Genomics', 'Genomic', 'Microarray','SNP chip','VCF', 'Linkage disequilibrium', 'Hardy-Weinberg equilibrium'],
  install_requires=['scipy'],
  classifiers=[
    'Development Status :: 4 - Beta',      #"3 - Alpha", "4 - Beta" or "5 - Production/Stable"
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3.0',   # Again, pick a license
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.9',
  ],
)
