from setuptools import setup, find_packages
import versioneer


def read_version():
    version = "none-0.0.0"
    with open("molpro/_version.py", "r") as f:
        for line in f.readlines():
            if "__version__" in line:
                version = line.split("=")[-1].strip().strip(' "').strip("'")
    return version


setup(
    name="molpro",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef>=1.23.0", "numpy>=1.12", "regex"],
)
