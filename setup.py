from setuptools import setup, find_packages
import versioneer

setup(
    name="pymolpro",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef>=1.30.0", "numpy>=1.12", "regex", "scipy>=1.9"],
)
