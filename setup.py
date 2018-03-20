from setuptools import setup, find_packages

setup(
    name='mw_plot',
    version='0.1.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic:: Scientific / Engineering:: Astronomy'],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'mw_plot': ['*.png']},
    python_requires='>=3.6',
    install_requires=[
        'numpy', 'astropy', 'matplotlib'],
    url='https://github.com/henrysky/milkyway_plot',
    project_urls={
        "Bug Tracker": "https://github.com/henrysky/milkyway_plot/issues",
        "Documentation": "https://github.com/henrysky/milkyway_plot",
        "Source Code": "https://github.com/henrysky/milkyway_plot/",
    },
    license='MIT',
    author='Henry Leung',
    author_email='henrysky.leung@mail.utoronto.ca',
    description='A handy python script to plot things on a face-on milkyway',
)
