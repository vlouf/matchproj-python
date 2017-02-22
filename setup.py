from setuptools import setup, find_packages

try:
    from numpy.distutils.core import setup
except:
    from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('msgr')
    config.add_subpackage('msgr.core')
    config.add_subpackage('msgr.core.instruments')
    config.add_subpackage('msgr.core.io')

    return config

setup(
    name='msgr',
    version='0.5',
    description='Matching Satellite and Ground Radar',
    long_description=readme,
    author='Valentin Louf',
    author_email='valentin.louf@bom.gov.au',
    url='https://gitlab.bom.gov.au/vlouf/Matchproj-python',
    license=license,
    packages=find_packages(exclude=('config', 'docs')),
    configuration=configuration,
    scripts=['scripts/matchvol', 'scripts/generate_config_matchvol']
)
