from setuptools import find_packages, setup

setup(
    name='mojitoo',
    version='0.1.0',
    description='single cell multimodal integration',
    url='https://github.com/CostaLab/MOJITOO/pymojitoo',
    author='Mingbo Cheng',
    author_email='chengmingbo@gmail.com',
    license='BSD 2-clause',
    install_requires=['numpy',
                      'sklearn',
                      'scipy',
                      'statsmodels',
                      'scanpy',
                      'anndata',
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    packages=find_packages()
)

