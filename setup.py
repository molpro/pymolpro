from setuptools import setup, find_packages


def read_version():
    version = "none-0.0.0"
    with open("pysjef_molpro/_version.py", "r") as f:
        for line in f.readlines():
            if "__version__" in line:
                version = line.split("=")[-1].strip().strip(' "').strip("'")
    return version


setup(
    name="pysjef_molpro",
    version=read_version(),
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef>=0.0.1", "numpy>=1.12", "regex"],
)
