import setuptools

setuptools.setup(
    name='probabilistic_epidemic_forecasting',
    version='0.0.1',
    author='Olivier Lapointe-Gagné and Bastian Raulier',
    author_email='olivier.lapointe-gagne.1@ulaval.ca',
    description='Epidemic parameters inference project',
    include_package_data=True,
    packages=setuptools.find_packages(),
    package_data={'probabilistic_epidemic_forecasting': ['epidemic_data/*']}
)
