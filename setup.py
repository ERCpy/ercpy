#!python
from distutils.core import setup

setup(name='test_ercpy',
      version='1.0',
      author="Martial Duchamp" "Jan Caron",
      author_email="martial.duchampl@gmail.com" "j.caron@fz-juelich.de",
      py_modules=['ercpy.holography',
                  'ercpy.eelsedx',
                  'ercpy.utils'],
      long_description="""Really long text here."""
      )
