from setuptools import setup, find_packages
import versioneer


def read_version():
    version = "none-0.0.0"
    with open("pymolpro/_version.py", "r") as f:
        for line in f.readlines():
            if "__version__" in line:
                version = line.split("=")[-1].strip().strip(' "').strip("'")
    return version


setup(
    name="pymolpro",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef>=1.25.1", "numpy>=1.12", "regex"],
)
