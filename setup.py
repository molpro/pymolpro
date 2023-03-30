from setuptools import setup, find_packages
import versioneer

setup(
    name="pymolpro",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef>=1.32.0", "numpy>=1.12", "regex", "scipy>=1.9"],
    package_dir={'pymolpro': 'pymolpro'},
    package_data={'pymolpro': ['share/database/*.json']},
)
