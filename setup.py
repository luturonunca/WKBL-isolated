from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='wkbl',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="something",
    license="MIT",
    author="ANC",
    author_email='ksnfkanf',
    url='https://github.com/as/wkbl',
    packages=['wkbl'],
    entry_points={
        'console_scripts': [
            'wkbl=wkbl.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='wkbl',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
