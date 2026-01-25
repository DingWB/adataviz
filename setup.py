from pathlib import Path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

this_directory = Path(__file__).parent

# release new version
setup(
    name="adataviz",
	use_scm_version={'version_scheme': 'post-release',"local_scheme": "no-local-version"},
	setup_requires=['setuptools_scm'],
	description="A python package to visualize adata",
    author="Wubin Ding",
    author_email="ding.wu.bin.gm@gmail.com",
    url="https://github.com/DingWB/adataviz",
    packages=["adataviz"],  # src
    install_requires=["matplotlib","numpy","pandas>=1.3.5", 
                      "scipy","loguru","fire","anndata","scanpy",
                      "PyComplexHeatmap","seaborn","adjustText",
                      "plotly","ipywidgets"],
    include_package_data=True,

    entry_points={
            'console_scripts': [
                'adataviz=adataviz:main',
                ],
        }
)





