from setuptools import setup, find_packages
import os

with open(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.rst"),
    encoding="utf-8",
) as f:
    long_description = f.read()

setup(
    name="mw_plot",
    version="0.11.0",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    packages=find_packages(),
    include_package_data=True,
    package_data={"mw_plot": ["*.png", "*.jpg"]},
    python_requires=">=3.8",
    install_requires=["requests", "numpy", "astropy", "matplotlib", "Pillow"],
    url="https://github.com/henrysky/milkyway_plot",
    project_urls={
        "Bug Tracker": "https://github.com/henrysky/milkyway_plot/issues",
        "Documentation": "https://milkyway-plot.readthedocs.io/",
        "Source Code": "https://github.com/henrysky/milkyway_plot/",
    },
    license="MIT",
    author="Henry Leung",
    author_email="henrysky.leung@utoronto.ca",
    description="A handy python package to do plotting on a face-on/edge-on/allsky map milkyway with matplotlib and bokeh",
    long_description=long_description,
)
