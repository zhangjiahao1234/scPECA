from setuptools import setup, find_packages

setup(
    name='scPECA',
    version='1.0',
    author='Jiahao Zhang',
    author_email='zhangjiahao@amss.ac.cn',
    url='https://zhuanlan.zhihu.com/p/26159930',
    description=u'吃枣药丸',
    packages=['jujube_pill'],
    install_requires=[],
    entry_points={
        'console_scripts': [
            'jujube=jujube_pill:jujube',
            'pill=jujube_pill:pill'
        ]
    }
)