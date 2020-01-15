from setuptools import setup, find_packages


setup(
    name="pysjef_molpro",
    packages=find_packages(),
    license="MIT",
    install_requires=["pysjef", "numpy>=1.12", "regex"],
)
