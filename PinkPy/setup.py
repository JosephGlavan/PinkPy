from distutils.core import setup, Extension

extension_mod = Extension("pinkpy", ["pinkpymodule.cpp"], include_dirs=['/usr/lib/python2.7/dist-packages/numpy/core/include/numpy/'])
setup(name = "pinkpy", version='1.0', description='Generate 3D 1/f pink noise', author='Joseph J Glavan', author_email='j.glavan4@gmail.com', ext_modules = [extension_mod])
