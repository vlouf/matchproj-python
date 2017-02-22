from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='msgr',
    version='0.5',
    description='Matching Satellite and Ground Radar',
    long_description=readme,
    author='Valentin Louf',
    author_email='valentin.louf@bom.gov.au',
    url='https://gitlab.bom.gov.au/vlouf/Matchproj-python',
    license=license,
    packages=find_packages(exclude=('config', 'docs'))
    scripts=['msgr/matchvol', 'msgr/generate_config_file']
)
