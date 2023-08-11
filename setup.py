import setuptools

setuptools.setup(
    name='repgenr',
    version='1.0',
    author='jaclew',
    description='no_description',
    packages=['repgenr'],
    scripts=['repgenr/repgenr','repgenr/metadata.py','repgenr/genome.py','repgenr/glance.py','repgenr/derep.py','repgenr/derep_worker.py','repgenr/derep_stocker.py','repgenr/derep_summarize.py','repgenr/phylo.py','repgenr/tree2tax.py','repgenr/x2fa.py']
)
