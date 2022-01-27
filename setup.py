from distutils.core import setup
setup(
  name = 'scarv',         # How you named your package folder (MyLib)
  packages = ['scarv'],   # Chose the same as "name"
  version = '0.4',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Selective Constraint Against Rare Variants',   # Give a short description about your library
  author = 'Joep ten Wolde',                   # Type in your name
  author_email = 'j.tenwolde@icloud.com',      # Type in your E-Mail
  url = 'https://github.com/jtenwolde/SCARV',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/jtenwolde/SCARV/archive/refs/tags/v_04.tar.gz',    # I explain this later on
  keywords = ['genetics', 'selection'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'numpy',
          'pandas',
          'pyranges',
          'xgboost',
          'sklearn',
          'itertools',
          'keras',
          'keras_genomics',
          'pysam',
          'gzip',
          'pyBigWig'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
