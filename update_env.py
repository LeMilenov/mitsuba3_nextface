import os
import sys
import shutil

PATH = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(PATH, 'out/build/x64-Release/RelWithDebInfo')

PYPATH = os.path.abspath(__file__)
PYPATH = os.path.join(os.path.dirname(PYPATH), '../.facenext/Scripts')

print(f'Copying {PATH} to {PYPATH}')
shutil.copytree(PATH, PYPATH, dirs_exist_ok=True)

# modules
modules = ["drjit", "mitsuba"]
PYPATH = os.path.abspath(__file__)
PYPATH = os.path.join(os.path.dirname(PYPATH), '../.facenext/Lib/site-packages')
PATH = os.path.join(PATH, 'python')
for module in modules:
    print(f'Copying {PATH}/{module} to {PYPATH}/{module}')
    shutil.copytree(f'{PATH}/{module}', f'{PYPATH}/{module}', dirs_exist_ok=True)